"""
=======================================================================
  DIAGRAMME D'INTERACTION N-M  —  Eurocode 2 (EN 1992-1-1)
  Section RECTANGULAIRE ou CIRCULAIRE en béton armé
=======================================================================

Référence théorique : cours ENPC BAEP1 – Séance 4 (J.M. Jaeger, Setec TPI)

Méthode de construction (§13 du cours) :
  Le plan de déformation pivote successivement autour des 3 pivots EC2 :

  Pivot A : εst = εud = 45‰  (déformation limite acier tendu)
            Le béton peut atteindre au max εcu2
            → domaine traction / flexion simple fortement sous-armée

  Pivot B : εc,sup = εcu2 = 3.5‰  (déformation ultime béton comprimé)
            L'acier inférieur varie de +εud → 0
            → domaine flexion simple standard jusqu'à pivot B+C

  Pivot C : εc(y) = εc2 = 2.0‰  (compression uniforme)
            εc,sup = 3.5‰ → εc,inf = 0 (transition B→C)
            puis εc,sup = εc,inf = 2.0‰ (compression pure)
            → domaine flexion composée avec forte compression

  Compression pure : εc = εc2 partout → N_Rd,max
  Traction pure    : εs = εud partout → N_Rd,min

Lois de comportement :
  - Béton  : parabole-rectangle EC2 §3.1.7
  - Acier  : bilinéaire avec palier horizontal EC2 §3.2.7
  - γc = 1.5 / γs = 1.15

Sections disponibles :
  - Rectangulaire b×h, armatures symétriques ou non (As_inf, As_sup)
  - Circulaire Ø D, n barres uniformément réparties

Sorties :
  - Tracé N-M  [kN / kN·m]  (N en abscisse horizontale, M en ordonnée verticale)
  - Vérification des sollicitations de calcul
  - Impression des points caractéristiques + état de déformation pivot B

=======================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path   import Path
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import warnings

# Mémorise le dernier calcul effectué via interface_graphique(), pour que
# les fonctions d'export (tracer_section, tracer, tracer_interactif...)
# puissent être rappelées ensuite sans avoir à ressaisir tous les
# arguments (section_type, section_params, mat, arma, Ac, ...).
_ETAT_GUI = {}

# ═══════════════════════════════════════════════════════════════
# 1.  PARAMÈTRES MATÉRIAUX EC2
# ═══════════════════════════════════════════════════════════════

def get_material_params(fck_MPa, fyk_MPa=500.0, gamma_c=1.5, gamma_s=1.15):
    """
    Paramètres matériaux selon EC2 tableau 3.1 / §3.2.

    gamma_c, gamma_s : coefficients partiels de sécurité matériaux
      (valeurs usuelles EC2 : 1,5 pour le béton, 1,15 pour l'acier —
      modifiables si besoin, p.ex. combinaisons accidentelles).

    Retourne un dict avec toutes les grandeurs utiles.
    La valeur εud du cours ENPC est 45‰ (pivot A) ; on conserve aussi
    εud_acier = 22.7‰ pour la capacité de déformation de l'acier seul
    (branche inclinée non utilisée ici, on prend le palier horizontal).
    """
    gc = gamma_c ;  gs = gamma_s
    fcd = fck_MPa / gc
    fyd = fyk_MPa / gs
    Es  = 200_000.0          # MPa

    if fck_MPa <= 50:
        n       = 2.0
        eps_c2  = 2.0e-3
        eps_cu2 = 3.5e-3
    else:                    # hautes résistances EC2 éq. 3.14-3.16
        n       = 1.4 + 23.4 * ((90 - fck_MPa) / 100) ** 4
        eps_c2  = (2.0 + 0.085 * (fck_MPa - 50) ** 0.53) * 1e-3
        eps_cu2 = (2.6 + 35.0 * ((90 - fck_MPa) / 100) ** 4) * 1e-3

    eps_yd  = fyd / Es        # déformation de plastification acier
    eps_ud  = 45.0e-3         # déformation limite acier PIVOT A  (cours ENPC §13)

    return dict(fck=fck_MPa, fcd=fcd, fyk=fyk_MPa, fyd=fyd, Es=Es,
                n=n, eps_c2=eps_c2, eps_cu2=eps_cu2,
                eps_yd=eps_yd, eps_ud=eps_ud,
                gamma_c=gc, gamma_s=gs)


# ═══════════════════════════════════════════════════════════════
# 2.  LOIS DE COMPORTEMENT
# ═══════════════════════════════════════════════════════════════

def sigma_c(eps, mat):
    """
    Contrainte béton [MPa] — loi parabole-rectangle EC2 §3.1.7.
    Convention : eps > 0 = COMPRESSION.
    """
    fcd, n, ec2, ecu = mat["fcd"], mat["n"], mat["eps_c2"], mat["eps_cu2"]
    if   eps <= 0.0  : return 0.0
    elif eps <= ec2  : return fcd * (1.0 - (1.0 - eps / ec2) ** n)
    elif eps <= ecu  : return fcd
    else             : return 0.0          # au-delà de εcu2 : rupture


def sigma_s(eps, mat):
    """
    Contrainte acier [MPa] — loi bilinéaire EC2 §3.2.7.
    Convention algébrique : eps > 0 = COMPRESSION, eps < 0 = TRACTION.
    Palier plastique jusqu'à εud (= 45‰, pivot A).
    """
    Es, fyd, eud = mat["Es"], mat["fyd"], mat["eps_ud"]
    sg = np.sign(eps) if eps != 0 else 1.0
    ae = abs(eps)
    if   ae <= fyd / Es : return Es * eps     # domaine élastique
    elif ae <= eud      : return sg * fyd     # palier plastique
    else                : return 0.0          # hors pivot A (ne doit pas arriver)


# ═══════════════════════════════════════════════════════════════
# 3.  ÉTATS DE DÉFORMATION — PIVOTS A, B, C  (cours ENPC §13)
# ═══════════════════════════════════════════════════════════════

def build_deformation_states(mat, h, n_piv_A=80, n_piv_B=120, n_piv_C=40):
    """
    Génère la liste des états de déformation (eps_top, eps_bot) qui font
    tourner le plan de déformation autour des pivots A → B → C.

    Convention : eps_top = déformation fibre HAUTE (comprimée en général)
                 eps_bot = déformation fibre BASSE  (tendue en général)
                 eps > 0 = compression

    Pivot A (εst = εud = 45‰, eps_bot = -εud, εc,sup ∈ [0 → εcu2])
    Pivot B (εc,sup = εcu2, eps_bot ∈ [-εud → 0])
    Transition B→C (εc,sup = εcu2, εc,inf ∈ [0 → εcu2])
              → correspond au balayage B→C du cours ENPC
    Pivot C  (εc uniforme de εcu2 → εc2, section totalement comprimée)
    """
    ecu = mat["eps_cu2"]
    ec2 = mat["eps_c2"]
    eud = mat["eps_ud"]

    states = []

    # ── Pivot A : εbot = -εud (traction pure → pivot A+B) ─────────────────
    # εtop varie de 0 (traction pure, axe neutre à l'infini côté bas)
    # jusqu'à εcu2 (pivot A+B)
    for eps_top in np.linspace(0.0, ecu, n_piv_A, endpoint=False):
        states.append((eps_top, -eud))

    # ── Pivot B : εtop = εcu2, εbot ∈ [-εud → 0] ─────────────────────────
    for eps_bot in np.linspace(-eud, 0.0, n_piv_B, endpoint=False):
        states.append((ecu, eps_bot))

    # ── Transition B→C : εtop = εcu2, εbot ∈ [0 → εcu2] ─────────────────
    # (section de plus en plus comprimée, axe neutre sort par le bas)
    for eps_bot in np.linspace(0.0, ecu, n_piv_C, endpoint=False):
        states.append((ecu, eps_bot))

    # ── Pivot C : déformation uniforme de εcu2 → εc2 ─────────────────────
    # Le diagramme de déf. est uniforme (axe neutre hors section côté bas)
    # L'état pivot C pur = εc2 partout
    for eps_uni in np.linspace(ecu, ec2, n_piv_C, endpoint=True):
        states.append((eps_uni, eps_uni))

    return states


# ═══════════════════════════════════════════════════════════════
# 4.  INTÉGRATION SUR LA SECTION
# ═══════════════════════════════════════════════════════════════

def compute_NM(eps_top, eps_bot, H, fibres, armatures, mat):
    """
    Calcule (N [kN], M [kN·m]) par intégration numérique pour un état
    de déformation (eps_top, eps_bot) donné.

    eps_top : déformation fibre HAUTE  (+ = compression)
    eps_bot : déformation fibre BASSE  (+ = compression)
    H       : hauteur totale de la section [m]
    fibres  : liste de (y_i, b_i, dy_i)  — y depuis CDG, + = haut
    armatures: liste de (y_i, As_i)
    """
    def eps_at(y):
        # Interpolation linéaire : y = +H/2 → eps_top ; y = -H/2 → eps_bot
        return eps_top + (eps_bot - eps_top) * (H / 2.0 - y) / H

    N_c = M_c = 0.0
    for (y_i, b_i, dy_i) in fibres:
        e_i  = eps_at(y_i)
        s_i  = sigma_c(e_i, mat)           # MPa
        dA   = b_i * dy_i                  # m²
        N_c += s_i * dA
        M_c += s_i * dA * y_i

    N_s = M_s = 0.0
    for (y_i, As_i) in armatures:
        e_i  = eps_at(y_i)
        # On soustrait la contrainte béton déplacée au droit de la barre
        s_ci = sigma_c(max(e_i, 0.0), mat)
        s_si = sigma_s(e_i, mat)
        fs   = (s_si - s_ci) * As_i        # MPa·m²
        N_s += fs
        M_s += fs * y_i

    # MPa·m² × 1e3 = kN  (1 MPa·m² = 1 MN·m²/1000 … non)
    # 1 MPa = 1 N/mm² = 1e6 N/m²  ; MPa × m² = 1e6 N = 1e3 kN ✓
    N = (N_c + N_s) * 1e3           # kN
    M = (M_c + M_s) * 1e3           # kN·m
    return N, M


# ═══════════════════════════════════════════════════════════════
# 5.  DISCRÉTISATION DES SECTIONS
# ═══════════════════════════════════════════════════════════════

def fibres_rect(b, h, n_div=200):
    """Tranches horizontales pour section rectangulaire b×h."""
    dy = h / n_div
    return [(h / 2 - (k + 0.5) * dy,   b,  dy) for k in range(n_div)]


def fibres_circ(D, n_div=200):
    """Tranches horizontales pour section circulaire de diamètre D."""
    R  = D / 2.0
    dy = D / n_div
    fibres = []
    for k in range(n_div):
        y_i = R - (k + 0.5) * dy
        b_i = 2.0 * np.sqrt(max(0.0, R**2 - y_i**2))
        if b_i > 1e-9:
            fibres.append((y_i, b_i, dy))
    return fibres


def armatures_rect(h, c_inf, c_sup, As_inf, As_sup):
    """
    Section rectangulaire : deux nappes horizontales.
    y_inf = -h/2 + c_inf  (fibre basse, tendue en flexion positive)
    y_sup = +h/2 - c_sup  (fibre haute, comprimée)
    """
    y_inf = -h / 2.0 + c_inf
    y_sup = +h / 2.0 - c_sup
    arma = []
    if As_sup > 0: arma.append((y_sup, As_sup))
    if As_inf > 0: arma.append((y_inf, As_inf))
    return arma


def armatures_circ(D, c_enr, nb_barres, As_tot):
    """n barres uniformément réparties sur cercle de diam. D - 2·c."""
    Rs     = (D - 2.0 * c_enr) / 2.0
    As_bar = As_tot / nb_barres
    return [(Rs * np.sin(np.pi / 2 + 2 * np.pi * i / nb_barres), As_bar)
            for i in range(nb_barres)]


# ═══════════════════════════════════════════════════════════════
# 6.  CONSTRUCTION DU DIAGRAMME D'INTERACTION
# ═══════════════════════════════════════════════════════════════

def diagramme_interaction(section_type, section_params,
                          fck, fyk=500.0,
                          gamma_c=1.5, gamma_s=1.15,
                          n_div=200, n_piv_A=80, n_piv_B=150, n_piv_C=50):
    """
    Paramètres communs :
      section_type   : "rect" ou "circ"
      section_params :
        rect → dict(b, h, c_inf, c_sup, As_inf, As_sup)
        circ → dict(D, c_enr, nb_barres, As_tot)
      fck, fyk : résistances caractéristiques [MPa]
      gamma_c, gamma_s : coefficients partiels matériaux (défaut EC2 usuel :
        1,5 / 1,15)

    Retourne : N_arr, M_arr [kN, kN·m], mat, pts_cles, fibres, armatures
    """
    mat = get_material_params(fck, fyk, gamma_c, gamma_s)

    if section_type == "rect":
        b, h = section_params["b"], section_params["h"]
        fibres   = fibres_rect(b, h, n_div)
        arma     = armatures_rect(h,
                                  section_params["c_inf"],
                                  section_params["c_sup"],
                                  section_params["As_inf"],
                                  section_params["As_sup"])
        H = h
        Ac = b * h

    elif section_type == "circ":
        D = section_params["D"]
        fibres = fibres_circ(D, n_div)
        arma   = armatures_circ(D,
                                section_params["c_enr"],
                                section_params["nb_barres"],
                                section_params["As_tot"])
        H  = D
        Ac = np.pi * (D / 2) ** 2
    else:
        raise ValueError("section_type doit être 'rect' ou 'circ'")

    # États de déformation (demi-enveloppe M ≥ 0)
    states = build_deformation_states(mat, H, n_piv_A, n_piv_B, n_piv_C)

    N_half, M_half = [], []

    # Traction pure (point de départ bas du diagramme)
    # NB : à εs = εud = 45‰ l'acier est très au-delà de εyd, donc σs = fyd
    # exactement (palier plastique) -> ce point est rigoureusement (N0, M0).
    eud = mat["eps_ud"]
    N0, M0 = compute_NM(-eud, -eud, H, fibres, arma, mat)
    N_half.append(N0) ;  M_half.append(M0)     # (ne plus forcer M=0)

    for (et, eb) in states:
        N, M = compute_NM(et, eb, H, fibres, arma, mat)
        N_half.append(N) ;  M_half.append(M)

    # Compression pure (point haut) : déformation uniforme εc2 = 2‰.
    # ATTENTION (cf. cours ENPC BAEP1 p.27-28, tableau "Poteau 60*60 C50") :
    # à εc2 = 2‰, l'acier n'est PAS plastifié pour du S500 (εyd = fyd/Es ≈
    # 2,174‰ > εc2) : le cours calcule explicitement σsc = Es·εc2 = 400 MPa
    # et NON σsc = fyd = 434,8 MPa. compute_NM() applique déjà correctement
    # sigma_s(eps, mat), qui redonne bien ce comportement élastique -> on
    # utilise directement (Np, Mp) sans passer par une formule fyd séparée.
    ec2 = mat["eps_c2"]
    Np, Mp = compute_NM(ec2, ec2, H, fibres, arma, mat)
    N_half.append(Np) ;  M_half.append(Mp)     # (ne plus forcer M=0)

    # ── Second demi-diagramme (l'"autre côté" M) ────────────────────────────
    # Construction rigoureuse pour un ferraillage QUELCONQUE (symétrique ou
    # non) : on réévalue les MÊMES états de déformation en échangeant les
    # rôles haut/bas (et, eb) -> (eb, et), SUR LA GÉOMÉTRIE RÉELLE FIXE
    # (fibres et armatures à leurs positions y réelles, non permutées).
    #
    # Justification : eps_at(y) est linéaire entre eps_top (en y=+H/2) et
    # eps_bot (en y=-H/2). En échangeant (et,eb), on montre que le nouveau
    # champ vérifie eps_at'(y) = eps_at(-y) : c'est exactement le champ de
    # déformation reflété par rapport au centre de la section. On intègre
    # ensuite ce champ reflété sur les positions RÉELLES (non reflétées)
    # des fibres béton et des armatures, ce qui donne le véritable second
    # état d'équilibre de la section fixe (chargement fléchissant dans
    # l'autre sens), sans jamais supposer As_sup = As_inf.
    #
    # NB : l'ancienne version fermait le diagramme par une symétrie miroir
    # M -> -M des mêmes points déjà calculés, ce qui n'est exact QUE si la
    # section est elle-même symétrique (As_sup=As_inf, c_sup=c_inf) — dans
    # le cas contraire cela produisait un contour incorrect (auto-croisé),
    # comme observé sur un ferraillage à nappe unique. La formule ci-dessous
    # se réduit exactement à l'ancienne symétrie quand la section EST
    # symétrique (elle ne change donc rien aux résultats déjà validés dans
    # le rapport de benchmark), et devient correcte dans le cas général.
    N_half2, M_half2 = [], []
    for (et, eb) in states:
        N, M = compute_NM(eb, et, H, fibres, arma, mat)
        N_half2.append(N) ;  M_half2.append(M)

    N_arr = np.array([N0] + N_half[1:-1] + [Np] + list(reversed(N_half2)) + [N0])
    M_arr = np.array([M0] + M_half[1:-1] + [Mp] + list(reversed(M_half2)) + [M0])

    # ── Repérage des zones de pivot (pour coloriage A/B/C du tracé) ─────────
    # Suit exactement l'ordre de construction ci-dessus : la zone "C" réunit
    # ici la transition B→C et le balayage pivot C pur (le cours et la
    # feuille H.Thonier ne distinguent que 3 couleurs sur le tracé).
    zones_states = (["A"] * n_piv_A + ["B"] * n_piv_B
                     + ["C"] * n_piv_C + ["C"] * n_piv_C)
    zones = (["A"] + zones_states + ["C"]
             + list(reversed(zones_states)) + ["A"])

    # ── Points caractéristiques ────────────────────────────────────────────
    # Nmax / Nmin proviennent directement de l'intégration fibre (Np, N0
    # ci-dessus) et non plus d'une formule séparée supposant à tort σsc=fyd
    # en compression pure (cf. cours ENPC p.27-28).
    As_tot = sum(a[1] for a in arma)
    Nmax   = Np
    Nmin   = N0
    idx_Mmax = int(np.argmax(M_arr))
    idx_Mmin = int(np.argmin(M_arr))
    Mmax     = float(M_arr[idx_Mmax])
    Mmin     = float(M_arr[idx_Mmin])
    N_atMmax = float(N_arr[idx_Mmax])
    N_atMmin = float(N_arr[idx_Mmin])

    # Point pivot B pur (εcu2 en haut, εyd en bas) — un des deux côtés ;
    # pour une section dissymétrique, l'autre côté (εyd en haut, εcu2 en
    # bas) peut différer et n'est pas reporté séparément ici.
    eps_yd = mat["eps_yd"]
    ecu    = mat["eps_cu2"]
    N_B, M_B = compute_NM(ecu, -eps_yd, H, fibres, arma, mat)

    pts_cles = {
        "N_Rd,max  (comp. pure)   [kN]"  : Nmax,
        "M à N_Rd,max             [kN·m]": Mp,
        "N_Rd,min  (tract. pure)  [kN]"  : Nmin,
        "M à N_Rd,min             [kN·m]": M0,
        "M_Rd,max                 [kN·m]": Mmax,
        "N au M_Rd,max            [kN]"  : N_atMmax,
        "M_Rd,min                 [kN·m]": Mmin,
        "N au M_Rd,min            [kN]"  : N_atMmin,
        "N au pivot B pur         [kN]"  : N_B,
        "M au pivot B pur         [kN·m]": M_B,
        "Aire section             [cm²]" : Ac * 1e4,
        "As totale                [cm²]" : As_tot * 1e4,
        "Taux d'armature ρ        [%]"   : As_tot / Ac * 100,
    }

    return N_arr, M_arr, mat, pts_cles, fibres, arma, Ac, H, zones


# ═══════════════════════════════════════════════════════════════
# 7.  TRACÉ
# ═══════════════════════════════════════════════════════════════

def dessiner_section(ax_sec, ax_info, section_type, section_params, mat, arma, Ac):
    """
    Dessine le schéma de section (rectangulaire ou circulaire, avec
    ferraillage) sur ax_sec, et les 3 encarts récap (Béton/Acier/
    Ferraillage) sur ax_info. Fonction partagée par tracer() (export PNG
    complet) et tracer_section() (export de la seule coupe).
    """
    As_tot = sum(a[1] for a in arma)

    ax_sec.set_xlim(-0.65, 0.65)
    ax_sec.set_ylim(-0.72, 0.72)
    ax_sec.set_aspect("equal")
    ax_sec.axis("off")
    ax_sec.set_title("Section", fontsize=9, fontweight="bold", pad=4)

    def _card(ax, x, y_top, title, body, edge_color, face_color):
        """Petit encart titré, bien délimité (bordure colorée + séparateur)."""
        txt = f"{title}\n{'─'*16}\n{body}"
        ax.text(x, y_top, txt, fontsize=6.6, va="top", family="monospace",
                color="#222", linespacing=1.4,
                bbox=dict(boxstyle="round,pad=0.30", fc=face_color,
                          ec=edge_color, lw=1.2, alpha=0.97))

    if section_type == "rect":
        b_s, h_s = section_params["b"], section_params["h"]
        c_inf_s, c_sup_s = section_params["c_inf"], section_params["c_sup"]
        As_inf_s, As_sup_s = section_params["As_inf"], section_params["As_sup"]
        sc = 0.46 / max(b_s, h_s)
        bS, hS = b_s * sc, h_s * sc
        rect = mpatches.FancyBboxPatch((-bS/2, -hS/2), bS, hS,
                                        boxstyle="square,pad=0",
                                        lw=1.5, ec="#333", fc="#cdd9ee")
        ax_sec.add_patch(rect)

        # Armatures — les DEUX nappes sont toujours représentées à leur
        # position réelle (y_sup, y_inf), même quand une nappe est vide :
        # dans ce cas on trace un repère en tirets plutôt que de simplement
        # omettre la nappe (ce qui pouvait laisser croire à un oubli plutôt
        # qu'à un choix de ferraillage volontaire).
        r_b = 0.022
        y_sup = (+h_s/2 - c_sup_s) * sc
        y_inf = (-h_s/2 + c_inf_s) * sc

        def _draw_nappe(y_pos, As_val, nb_reel):
            if As_val > 1e-9:
                # Nombre de points affichés = nombre RÉEL de barres saisi
                # dans le formulaire (nb_sup/nb_inf), s'il est disponible.
                # À défaut (appel direct de dessiner_section sans cette
                # info), on retombe sur une estimation à partir de la
                # surface — approximative, car elle suppose un diamètre
                # de barre de référence (16mm) qui peut différer du réel.
                if nb_reel is not None:
                    nb = max(1, int(nb_reel))
                else:
                    nb = max(2, min(8, round(As_val / (np.pi * (0.016/2)**2))))
                for xb in np.linspace(-bS/2 + 0.03, bS/2 - 0.03, nb):
                    ax_sec.add_patch(plt.Circle((xb, y_pos), r_b,
                                                color="#b71c1c", zorder=4))
            else:
                ax_sec.plot([-bS/2 + 0.03, bS/2 - 0.03], [y_pos, y_pos],
                            color="#b0b0b0", lw=1.2, ls=(0, (4, 3)), zorder=4)

        _draw_nappe(y_sup, As_sup_s, section_params.get("nb_sup"))
        _draw_nappe(y_inf, As_inf_s, section_params.get("nb_inf"))

        # Cotes
        ax_sec.annotate("", xy=(bS/2+0.07,-hS/2), xytext=(bS/2+0.07, hS/2),
                        arrowprops=dict(arrowstyle="<->",color="#333",lw=1.0))
        ax_sec.text(bS/2+0.14, 0, f"h={h_s*100:.0f}cm",
                    va="center", fontsize=7.5, rotation=90)
        ax_sec.annotate("", xy=(-bS/2,-hS/2-0.07), xytext=(bS/2,-hS/2-0.07),
                        arrowprops=dict(arrowstyle="<->",color="#333",lw=1.0))
        ax_sec.text(0, -hS/2-0.14, f"b={b_s*100:.0f}cm",
                    ha="center", fontsize=7.5)

        box_beton = (f" fck=C{mat['fck']:.0f}   fcd={mat['fcd']:.1f} MPa\n"
                     f" εcu2={mat['eps_cu2']*1e3:.1f}‰   εc2={mat['eps_c2']*1e3:.1f}‰")
        box_acier = (f" fyk=S{int(mat['fyk'])}   fyd={mat['fyd']:.1f} MPa\n"
                     f" εud={mat['eps_ud']*1e3:.1f}‰ (pivot A)")
        box_fer   = (f" As_sup = {As_sup_s*1e4:5.1f} cm²  (c={c_sup_s*100:.1f}cm)\n"
                     f" As_inf = {As_inf_s*1e4:5.1f} cm²  (c={c_inf_s*100:.1f}cm)\n"
                     f" ρ total = {As_tot/Ac*100:.2f} %")

    else:  # circulaire
        D_s  = section_params["D"]
        c_e  = section_params["c_enr"]
        sc   = 0.46 / (D_s / 2)
        R_d  = D_s / 2 * sc
        Rs_d = (D_s / 2 - c_e) * sc
        ax_sec.add_patch(plt.Circle((0,0), R_d,  fc="#cdd9ee", ec="#333", lw=1.5))
        ax_sec.add_patch(plt.Circle((0,0), Rs_d, fill=False, ec="#999",
                                     ls="--", lw=0.8))
        nb = section_params["nb_barres"]
        r_b = 0.024
        for i in range(nb):
            a = np.pi/2 + 2*np.pi*i/nb
            ax_sec.add_patch(plt.Circle((Rs_d*np.cos(a), Rs_d*np.sin(a)),
                                         r_b, color="#b71c1c", zorder=4))
        ax_sec.annotate("", xy=(R_d+0.07,-R_d), xytext=(R_d+0.07, R_d),
                        arrowprops=dict(arrowstyle="<->",color="#333",lw=1.0))
        ax_sec.text(R_d+0.14, 0, f"Ø{D_s*100:.0f}cm",
                    va="center", fontsize=7.5, rotation=90)

        box_beton = (f" fck=C{mat['fck']:.0f}   fcd={mat['fcd']:.1f} MPa\n"
                     f" εcu2={mat['eps_cu2']*1e3:.1f}‰   εc2={mat['eps_c2']*1e3:.1f}‰")
        box_acier = (f" fyk=S{int(mat['fyk'])}   fyd={mat['fyd']:.1f} MPa\n"
                     f" εud={mat['eps_ud']*1e3:.1f}‰ (pivot A)")
        box_fer   = (f" {nb} barres  (enrobage {c_e*100:.1f}cm)\n"
                     f" As_tot = {As_tot*1e4:.1f} cm²\n"
                     f" ρ total = {As_tot/Ac*100:.2f} %")

    # Récapitulatif en 3 encarts bien séparés
    _card(ax_info, 0.04, 0.97, "BÉTON",       box_beton, "#1565C0", "#E3F2FD")
    _card(ax_info, 0.04, 0.65, "ACIER",       box_acier, "#EF6C00", "#FFF3E0")
    _card(ax_info, 0.04, 0.35, "FERRAILLAGE", box_fer,   "#2E7D32", "#E8F5E9")


def tracer_section(section_type=None, section_params=None, mat=None,
                    arma=None, Ac=None, nom_fichier="section_EC2.png"):
    """
    Exporte UNIQUEMENT la coupe de section (schéma + encarts Béton/Acier/
    Ferraillage), sans le diagramme N-M — utile en complément du diagramme
    interactif Plotly (tracer_interactif), qui ne trace pas la section.

    Appelée sans argument, reprend automatiquement les valeurs du dernier
    calcul lancé depuis interface_graphique() (bouton "Calculer") — pratique
    pour ré-exporter la coupe dans une cellule séparée sans tout ressaisir :

        ec2.interface_graphique()   # ... on clique sur Calculer ...
        ec2.tracer_section()        # ré-exporte la coupe, mêmes données
    """
    if section_type is None:
        if not _ETAT_GUI:
            raise RuntimeError(
                "Aucun calcul disponible : lancez d'abord interface_graphique() "
                "et cliquez sur \"Calculer\", ou passez explicitement "
                "section_type, section_params, mat, arma, Ac.")
        section_type   = _ETAT_GUI["section_type"]
        section_params = _ETAT_GUI["section_params"]
        mat            = _ETAT_GUI["mat"]
        arma           = _ETAT_GUI["arma"]
        Ac             = _ETAT_GUI["Ac"]

    fig = plt.figure(figsize=(4.6, 8.6), facecolor="#f4f6f9")
    gs = GridSpec(2, 1, height_ratios=[1.15, 1.85], hspace=0.08, figure=fig)
    ax_sec  = fig.add_subplot(gs[0])
    ax_info = fig.add_subplot(gs[1])
    ax_info.set_xlim(0, 1) ;  ax_info.set_ylim(0, 1) ;  ax_info.axis("off")

    dessiner_section(ax_sec, ax_info, section_type, section_params, mat, arma, Ac)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="This figure includes Axes")
        plt.tight_layout()
    plt.savefig(nom_fichier, dpi=150, bbox_inches="tight")
    print(f"  Coupe de section sauvegardée : {nom_fichier}")
    plt.show()


# ═══════════════════════════════════════════════════════════════
# 7ter.  EXPORT CSV — coordonnées du diagramme d'interaction
# ═══════════════════════════════════════════════════════════════

def exporter_csv(N_arr=None, M_arr=None, zones=None,
                  nom_fichier="diagramme_interaction_EC2.csv"):
    """
    Exporte les coordonnées (N, M) du contour du diagramme d'interaction
    au format CSV — une ligne par point, dans l'ordre de tracé, avec la
    zone de pivot (A/B/C) associée en 3e colonne.

    Format : séparateur ';' et virgule décimale (convention Excel FR).

    Appelée sans argument, reprend automatiquement les valeurs du dernier
    calcul lancé depuis interface_graphique().
    """
    if N_arr is None:
        if not _ETAT_GUI:
            raise RuntimeError(
                "Aucun calcul disponible : lancez d'abord interface_graphique() "
                "et cliquez sur \"Calculer\", ou passez explicitement N_arr, M_arr.")
        N_arr = _ETAT_GUI["N_arr"] ;  M_arr = _ETAT_GUI["M_arr"]
        zones = _ETAT_GUI["zones"]

    if not nom_fichier.lower().endswith(".csv"):
        nom_fichier += ".csv"

    avec_zone = zones is not None and len(zones) == len(N_arr)

    import csv
    with open(nom_fichier, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f, delimiter=";")
        entete = ["N_Ed [kN]", "M_Ed [kN.m]"]
        if avec_zone:
            entete.append("Pivot")
        writer.writerow(entete)
        for i in range(len(N_arr)):
            ligne = [f"{N_arr[i]:.3f}".replace(".", ","),
                     f"{M_arr[i]:.3f}".replace(".", ",")]
            if avec_zone:
                ligne.append(zones[i])
            writer.writerow(ligne)

    print(f"  Coordonnées exportées : {nom_fichier}  ({len(N_arr)} points)")


# ═══════════════════════════════════════════════════════════════
# 7quater.  EXPORT PDF — rapport complet (hypothèses + diagramme)
# ═══════════════════════════════════════════════════════════════

def exporter_pdf(N_arr=None, M_arr=None, mat=None, pts_cles=None,
                  section_type=None, section_params=None,
                  fibres=None, arma=None, Ac=None, H=None,
                  zones=None, sollicitations=None,
                  nom_fichier="rapport_diagramme_interaction_EC2.pdf"):
    """
    Exporte le rapport complet au format PDF : hypothèses de calcul
    (section, matériaux béton/acier, ferraillage — les mêmes encarts que
    sur l'export PNG) accompagnées du diagramme d'interaction N-M.

    Reprend exactement la mise en page de tracer(), simplement enregistrée
    en PDF (vectoriel) plutôt qu'en PNG (matriciel) — utile pour joindre
    à une note de calcul.

    Appelée sans argument, reprend automatiquement les valeurs du dernier
    calcul lancé depuis interface_graphique().
    """
    if not nom_fichier.lower().endswith(".pdf"):
        nom_fichier += ".pdf"

    tracer(N_arr, M_arr, mat, pts_cles,
           section_type=section_type, section_params=section_params,
           fibres=fibres, arma=arma, Ac=Ac, H=H,
           zones=zones, sollicitations=sollicitations,
           nom_fichier=nom_fichier)


def tracer(N_arr=None, M_arr=None, mat=None, pts_cles=None,
           section_type=None, section_params=None,
           fibres=None, arma=None, Ac=None, H=None,
           zones=None,
           sollicitations=None,
           nom_fichier="diagramme_interaction_EC2.png"):
    """
    Diagramme d'interaction N-M complet (courbe + schéma de section +
    encarts Béton/Acier/Ferraillage), exporté au format déduit de
    l'extension de nom_fichier (.png, .pdf, .svg... — tout ce que
    matplotlib sait écrire).

    Appelée sans argument, reprend automatiquement les valeurs du dernier
    calcul lancé depuis interface_graphique().
    """
    if N_arr is None:
        if not _ETAT_GUI:
            raise RuntimeError(
                "Aucun calcul disponible : lancez d'abord interface_graphique() "
                "et cliquez sur \"Calculer\", ou passez explicitement tous les "
                "arguments.")
        N_arr = _ETAT_GUI["N_arr"] ;              M_arr = _ETAT_GUI["M_arr"]
        mat = _ETAT_GUI["mat"] ;                  pts_cles = _ETAT_GUI["pts_cles"]
        section_type = _ETAT_GUI["section_type"]
        section_params = _ETAT_GUI["section_params"]
        fibres = _ETAT_GUI["fibres"] ;            arma = _ETAT_GUI["arma"]
        Ac = _ETAT_GUI["Ac"] ;                    H = _ETAT_GUI["H"]
        zones = _ETAT_GUI["zones"]
        if sollicitations is None:
            sollicitations = _ETAT_GUI.get("sollicitations")

    fcd = mat["fcd"] ;  fyd = mat["fyd"]
    As_tot = sum(a[1] for a in arma)

    ncols = 2
    fig = plt.figure(figsize=(5 * ncols + 2, 10.5), facecolor="#f4f6f9")
    gs  = GridSpec(2, ncols,
                   width_ratios=[3.0, 1.2],
                   height_ratios=[5.6, 1.0],
                   wspace=0.14, hspace=0.28, figure=fig)

    ax_nm = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[1, :])
    ax_leg.axis("off")

    # Colonne de droite scindée en 2 zones INDÉPENDANTES : le dessin de la
    # section (petit, en haut) et les 3 encarts récap (grands, en bas).
    # Les avoir sur deux axes séparés (au lieu d'empiler du texte sur un
    # seul axe à la main) évite tout chevauchement, quel que soit le
    # nombre de lignes de chaque encart.
    gs_right = gs[0, 1].subgridspec(2, 1, height_ratios=[1.15, 1.85], hspace=0.08)
    ax_sec  = fig.add_subplot(gs_right[0])
    ax_info = fig.add_subplot(gs_right[1])
    ax_info.set_xlim(0, 1) ;  ax_info.set_ylim(0, 1) ;  ax_info.axis("off")

    palette_sol = ["#E65100","#1B5E20","#4A148C","#006064","#BF360C","#880E4F"]

    # Couleurs par zone de pivot (convention du cours / feuille H.Thonier)
    zone_colors = {"A": "#C62828",   # rouge  — pivot A (traction)
                   "B": "#1565C0",   # bleu   — pivot B (flexion)
                   "C": "#2E7D32"}   # vert   — pivot C (compression)
    zone_labels = {"A": "Pivot A (traction, εs=εud=45‰)",
                   "B": "Pivot B (flexion, εc,sup=εcu2=3,5‰)",
                   "C": "Pivot C (compression, transition B→C + εc2=2‰)"}

    # ─── Diagramme N-M ──────────────────────────────────────────────────────
    def _plot_diag(ax, xs, ys, zns, xlabel, ylabel, title_suffix,
                  sols_x, sols_y, sol_labels):
        ax.fill(xs, ys, alpha=0.10, color="#1565C0", zorder=1)

        if zns is not None and len(zns) == len(xs):
            pts  = np.column_stack([xs, ys])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            cols = [zone_colors.get(z, "#1565C0") for z in zns[:-1]]
            lc = LineCollection(segs, colors=cols, linewidths=2.2, zorder=3)
            ax.add_collection(lc)
            handles = [Line2D([0], [0], color=zone_colors[z], lw=2.2,
                               label=zone_labels[z]) for z in ("A", "B", "C")]
        else:
            ax.plot(xs, ys, color="#1565C0", lw=2.2, zorder=3,
                    label="Diagramme d'interaction EC2")
            handles = []

        ax.axhline(0, color="#000", lw=0.6, ls="--", alpha=0.4)
        ax.axvline(0, color="#000", lw=0.6, ls="--", alpha=0.4)

        # Point M_Rd,max (M porté en ordonnée ys)
        idx = int(np.argmax(ys))
        h_mmax, = ax.plot(xs[idx], ys[idx], "o", color="#B71C1C", ms=7, zorder=5,
                label=f"M_Rd,max = {ys[idx]:.1f}")
        handles.append(h_mmax)

        # Pivot B pur : N en abscisse, M en ordonnée
        N_piv_B = pts_cles["N au pivot B pur         [kN]"]
        M_piv_B = pts_cles["M au pivot B pur         [kN·m]"]
        h_pivB, = ax.plot(N_piv_B, M_piv_B, "s", color="#FF8F00", ms=7, zorder=5,
                label=f"Pivot B pur  N={N_piv_B:.1f} / M={M_piv_B:.1f}")
        handles.append(h_pivB)

        # Sollicitations
        if sols_x:
            for i, (sx, sy, lbl) in enumerate(zip(sols_x, sols_y, sol_labels)):
                c = palette_sol[i % len(palette_sol)]
                h, = ax.plot(sx, sy, "D", color=c, ms=9, zorder=6,
                        label=f"{lbl} ({sx:.1f} ; {sy:.1f})")
                handles.append(h)
                ax.annotate(lbl, (sx, sy), xytext=(7, 7),
                            textcoords="offset points", fontsize=8, color=c,
                            fontweight="bold")

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title_suffix, fontsize=10, fontweight="bold")
        ax.grid(True, alpha=0.3)
        ax.set_facecolor("#fdfeff")
        ax.relim(); ax.autoscale_view()
        return handles

    # Axes N-M
    titre_section = (
        f"Section RECTANGULAIRE {section_params['b']*100:.0f}×"
        f"{section_params['h']*100:.0f} cm"
        if section_type == "rect"
        else f"Section CIRCULAIRE Ø{section_params['D']*100:.0f} cm"
    )
    titre = (f"Diagramme N-M — {titre_section}\n"
             f"C{mat['fck']:.0f} / S{mat['fyk']:.0f} — "
             f"As={As_tot*1e4:.1f} cm²  ρ={As_tot/Ac*100:.2f}%")

    sols_M = [s[1] for s in sollicitations] if sollicitations else []
    sols_N = [s[0] for s in sollicitations] if sollicitations else []
    sols_L = [s[2] for s in sollicitations] if sollicitations else []

    # N en abscisse (horizontal), M en ordonnée (vertical)
    handles = _plot_diag(ax_nm, N_arr, M_arr, zones,
               "N$_{Ed}$ [kN]", "M$_{Ed}$ [kN·m]", titre,
               sols_N, sols_M, sols_L)

    ax_leg.legend(handles=handles, loc="center", ncol=3, fontsize=7.8,
                  frameon=True, framealpha=0.95, edgecolor="#999",
                  borderaxespad=0.2, handletextpad=0.6, columnspacing=1.4)

    # ─── Schéma de la section ───────────────────────────────────────────────
    dessiner_section(ax_sec, ax_info, section_type, section_params, mat, arma, Ac)

    # ─── Console ────────────────────────────────────────────────────────────
    print("\n" + "═"*60)
    print("  POINTS CARACTÉRISTIQUES  —  EC2 / ENPC BAEP1")
    print("═"*60)
    for k, v in pts_cles.items():
        print(f"  {k:<42} {v:>10.1f}")
    print("═"*60)

    if sollicitations:
        print("\n  VÉRIFICATION DES SOLLICITATIONS DE CALCUL")
        print("  " + "─"*56)
        path = Path(list(zip(M_arr, N_arr)))
        for (Ned, Med, lbl) in sollicitations:
            ok  = path.contains_point((Med, Ned))
            tag = "✔  VÉRIFIÉ    " if ok else "✘  DÉPASSEMENT"
            print(f"  {lbl:14s}  NEd={Ned:8.0f}kN  MEd={Med:7.0f}kN·m  →  {tag}")
        print("  " + "─"*56)

    plt.suptitle("Diagramme d'interaction  —  Eurocode 2  (pivots A / B / C)",
                 fontsize=11, fontweight="bold", y=1.01)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="This figure includes Axes")
        plt.tight_layout()
    plt.savefig(nom_fichier, dpi=150, bbox_inches="tight")
    print(f"\n  Figure sauvegardée : {nom_fichier}")
    plt.show()


# ═══════════════════════════════════════════════════════════════
# 7bis.  VERSION INTERACTIVE (Plotly) — coordonnées au survol
# ═══════════════════════════════════════════════════════════════

def tracer_interactif(N_arr=None, M_arr=None, mat=None, pts_cles=None,
                       section_type=None, section_params=None,
                       fibres=None, arma=None, Ac=None, H=None,
                       zones=None,
                       sollicitations=None,
                       nom_fichier="diagramme_interaction_EC2.html",
                       auto_open=False):
    """
    Version interactive du diagramme N-M (Plotly) : en promenant la souris
    sur la courbe ou sur un point, les coordonnées (N, M) s'affichent en
    infobulle. Produit un fichier .html autonome (ouvrable dans n'importe
    quel navigateur, sans connexion internet) et, dans un notebook Jupyter,
    s'affiche aussi directement dans la cellule.

    Reprend exactement les mêmes données que tracer() (mêmes N_arr, M_arr,
    zones, points caractéristiques) — seul le rendu change.

    Appelée sans argument, reprend automatiquement les valeurs du dernier
    calcul lancé depuis interface_graphique().
    """
    if N_arr is None:
        if not _ETAT_GUI:
            raise RuntimeError(
                "Aucun calcul disponible : lancez d'abord interface_graphique() "
                "et cliquez sur \"Calculer\", ou passez explicitement tous les "
                "arguments.")
        N_arr = _ETAT_GUI["N_arr"] ;              M_arr = _ETAT_GUI["M_arr"]
        mat = _ETAT_GUI["mat"] ;                  pts_cles = _ETAT_GUI["pts_cles"]
        section_type = _ETAT_GUI["section_type"]
        section_params = _ETAT_GUI["section_params"]
        fibres = _ETAT_GUI["fibres"] ;            arma = _ETAT_GUI["arma"]
        Ac = _ETAT_GUI["Ac"] ;                    H = _ETAT_GUI["H"]
        zones = _ETAT_GUI["zones"]
        if sollicitations is None:
            sollicitations = _ETAT_GUI.get("sollicitations")

    import plotly.graph_objects as go

    As_tot = sum(a[1] for a in arma)
    zone_colors = {"A": "#C62828", "B": "#1565C0", "C": "#2E7D32"}
    zone_names  = {"A": "Pivot A (traction, εs=εud=45‰)",
                   "B": "Pivot B (flexion, εc,sup=εcu2=3,5‰)",
                   "C": "Pivot C (compression, transition B→C + εc2=2‰)"}

    fig = go.Figure()

    # Contour rempli (une seule trace, sans survol dédié : sert de fond)
    fig.add_trace(go.Scatter(
        x=N_arr, y=M_arr, mode="lines",
        line=dict(color="rgba(0,0,0,0)"),
        fill="toself", fillcolor="rgba(21,101,192,0.10)",
        hoverinfo="skip", showlegend=False))

    # Contour colorié par zone de pivot, AVEC coordonnées au survol
    if zones is not None and len(zones) == len(N_arr):
        for z in ("A", "B", "C"):
            idx = [i for i, zz in enumerate(zones) if zz == z]
            # on découpe en segments contigus pour ne pas relier des tronçons
            # non adjacents par une ligne parasite
            segs, cur = [], []
            for i in range(len(zones)):
                if zones[i] == z:
                    cur.append(i)
                else:
                    if cur:
                        segs.append(cur); cur = []
            if cur:
                segs.append(cur)
            first = True
            for seg in segs:
                # inclut le point suivant pour la continuité visuelle du tracé
                idxs = seg + ([seg[-1] + 1] if seg[-1] + 1 < len(N_arr) else [])
                fig.add_trace(go.Scatter(
                    x=N_arr[idxs], y=M_arr[idxs], mode="lines",
                    line=dict(color=zone_colors[z], width=2.5),
                    name=zone_names[z], legendgroup=z,
                    showlegend=first,
                    hovertemplate="N = %{x:.1f} kN<br>M = %{y:.1f} kN·m<extra>"
                                  + zone_names[z] + "</extra>"))
                first = False
    else:
        fig.add_trace(go.Scatter(
            x=N_arr, y=M_arr, mode="lines",
            line=dict(color="#1565C0", width=2.5),
            name="Diagramme d'interaction EC2",
            hovertemplate="N = %{x:.1f} kN<br>M = %{y:.1f} kN·m<extra></extra>"))

    # Points caractéristiques
    N_piv_B = pts_cles["N au pivot B pur         [kN]"]
    M_piv_B = pts_cles["M au pivot B pur         [kN·m]"]
    idx_mmax = int(np.argmax(M_arr))
    fig.add_trace(go.Scatter(
        x=[N_arr[idx_mmax]], y=[M_arr[idx_mmax]], mode="markers",
        marker=dict(color="#B71C1C", size=11, symbol="circle"),
        name=f"M_Rd,max = {M_arr[idx_mmax]:.1f} kN·m",
        hovertemplate="M_Rd,max<br>N = %{x:.1f} kN<br>M = %{y:.1f} kN·m<extra></extra>"))
    fig.add_trace(go.Scatter(
        x=[N_piv_B], y=[M_piv_B], mode="markers",
        marker=dict(color="#FF8F00", size=11, symbol="square"),
        name=f"Pivot B pur (N={N_piv_B:.1f} / M={M_piv_B:.1f})",
        hovertemplate="Pivot B pur<br>N = %{x:.1f} kN<br>M = %{y:.1f} kN·m<extra></extra>"))

    # Sollicitations de calcul
    if sollicitations:
        path = Path(list(zip(M_arr, N_arr)))
        for (Ned, Med, lbl) in sollicitations:
            ok = path.contains_point((Med, Ned))
            fig.add_trace(go.Scatter(
                x=[Ned], y=[Med], mode="markers+text",
                marker=dict(color="#8E24AA" if ok else "#D81B60", size=13,
                            symbol="diamond",
                            line=dict(color="white", width=1)),
                text=[lbl], textposition="top center",
                name=f"{lbl} — {'vérifié' if ok else 'dépassement'}",
                hovertemplate=(f"{lbl}<br>N = %{{x:.1f}} kN<br>M = %{{y:.1f}} kN·m<br>"
                                f"{'✔ vérifié' if ok else '✘ dépassement'}<extra></extra>")))

    titre = (f"Diagramme N-M — Section "
             f"{'RECTANGULAIRE' if section_type=='rect' else 'CIRCULAIRE'} — "
             f"C{mat['fck']:.0f} / S{mat['fyk']:.0f} — "
             f"As={As_tot*1e4:.1f} cm²  ρ={As_tot/Ac*100:.2f}%")

    fig.update_layout(
        title=dict(text=titre, x=0.02, font=dict(size=15)),
        xaxis_title="N_Ed [kN]", yaxis_title="M_Ed [kN·m]",
        hovermode="closest",
        template="plotly_white",
        width=980, height=680,
        legend=dict(orientation="h", yanchor="bottom", y=-0.22,
                    xanchor="left", x=0),
    )
    fig.add_hline(y=0, line=dict(color="#999", width=1, dash="dot"))
    fig.add_vline(x=0, line=dict(color="#999", width=1, dash="dot"))

    fig.write_html(nom_fichier, auto_open=auto_open)
    print(f"  Figure interactive sauvegardée : {nom_fichier}")
    return fig


# ═══════════════════════════════════════════════════════════════
# 7ter.  INTERFACE GRAPHIQUE (ipywidgets) — pour notebook Jupyter
# ═══════════════════════════════════════════════════════════════

def interface_graphique():
    """
    Affiche un formulaire interactif (ipywidgets) dans une cellule Jupyter :
    type de section, matériaux (fck, fyk, γc, γs — préremplis aux valeurs
    usuelles EC2), géométrie et ferraillage. Un bouton "Calculer" lance
    diagramme_interaction() puis affiche le diagramme interactif Plotly.

    Nécessite le paquet ipywidgets : !pip install ipywidgets --quiet
    (puis redémarrer le kernel si c'est la première installation).
    """
    import ipywidgets as widgets
    from IPython.display import display, clear_output

    style = {"description_width": "140px"}
    layout_champ = widgets.Layout(width="280px")

    # ── Section & matériaux communs ─────────────────────────────────────
    w_forme = widgets.ToggleButtons(
        options=[("Rectangulaire", "rect"), ("Circulaire", "circ")],
        description="Type de section :", style=style)

    # Classes de béton normalisées — EC2 Tableau 3.1 (fck cylindre, MPa)
    CLASSES_BETON = [
        ("C12/15", 12.0), ("C16/20", 16.0), ("C20/25", 20.0),
        ("C25/30", 25.0), ("C30/37", 30.0), ("C35/45", 35.0),
        ("C40/50", 40.0), ("C45/55", 45.0), ("C50/60", 50.0),
        ("C55/67", 55.0), ("C60/75", 60.0), ("C70/85", 70.0),
        ("C80/95", 80.0), ("C90/105", 90.0),
    ]
    w_fck = widgets.Dropdown(options=CLASSES_BETON, value=25.0,
                              description="Classe béton :",
                              style=style, layout=layout_champ)
    w_fyk = widgets.FloatText(value=500.0, description="fyk [MPa] :",
                               style=style, layout=layout_champ)
    w_gc  = widgets.FloatText(value=1.5, description="γc :",
                               style=style, layout=layout_champ)
    w_gs  = widgets.FloatText(value=1.15, description="γs :",
                               style=style, layout=layout_champ)

    box_materiaux = widgets.VBox([
        widgets.HTML("<b>Matériaux</b>"),
        widgets.HBox([w_fck, w_fyk]),
        widgets.HBox([w_gc, w_gs]),
    ])

    # ── Rectangulaire ────────────────────────────────────────────────────
    w_b    = widgets.FloatText(value=1.00, description="b [m] :", style=style, layout=layout_champ)
    w_h    = widgets.FloatText(value=0.40, description="h [m] :", style=style, layout=layout_champ)
    w_cinf = widgets.FloatText(value=0.055, description="Enrobage inf. [m] :", style=style, layout=layout_champ)
    w_csup = widgets.FloatText(value=0.055, description="Enrobage sup. [m] :", style=style, layout=layout_champ)
    w_nbHAs = widgets.IntText(value=5,  description="Nb barres sup. :", style=style, layout=layout_champ)
    w_phis  = widgets.FloatText(value=10.0, description="Ø sup. [mm] :", style=style, layout=layout_champ)
    w_nbHAi = widgets.IntText(value=10, description="Nb barres inf. :", style=style, layout=layout_champ)
    w_phii  = widgets.FloatText(value=10.0, description="Ø inf. [mm] :", style=style, layout=layout_champ)

    box_rect = widgets.VBox([
        widgets.HTML("<b>Géométrie</b>"),
        widgets.HBox([w_b, w_h]),
        widgets.HBox([w_csup, w_cinf]),
        widgets.HTML("<b>Nappe supérieure</b>"),
        widgets.HBox([w_nbHAs, w_phis]),
        widgets.HTML("<b>Nappe inférieure</b>"),
        widgets.HBox([w_nbHAi, w_phii]),
    ])

    # ── Circulaire ───────────────────────────────────────────────────────
    w_D    = widgets.FloatText(value=0.60, description="D [m] :", style=style, layout=layout_champ)
    w_enr  = widgets.FloatText(value=0.070, description="Enrobage [m] :", style=style, layout=layout_champ)
    w_nbb  = widgets.IntText(value=12, description="Nb barres :", style=style, layout=layout_champ)
    w_phib = widgets.FloatText(value=16.0, description="Ø barres [mm] :", style=style, layout=layout_champ)

    box_circ = widgets.VBox([
        widgets.HTML("<b>Géométrie & ferraillage</b>"),
        widgets.HBox([w_D, w_enr]),
        widgets.HBox([w_nbb, w_phib]),
    ])

    box_section = widgets.VBox([box_rect])

    def _on_forme_change(change):
        box_section.children = [box_rect] if w_forme.value == "rect" else [box_circ]
    w_forme.observe(_on_forme_change, names="value")

    # ── Sollicitations à vérifier (facultatif) ──────────────────────────
    w_sol_txt = widgets.Textarea(
        value="", placeholder="Une ligne par cas, ex. : 2189, 842, ELU1",
        description="Sollicitations\n(N,M,label) :", style=style,
        layout=widgets.Layout(width="420px", height="70px"))

    btn = widgets.Button(description="Calculer", button_style="primary", icon="play")
    out = widgets.Output()

    def _on_click(_b):
        with out:
            clear_output(wait=True)
            try:
                forme = w_forme.value
                fck, fyk, gc, gs = w_fck.value, w_fyk.value, w_gc.value, w_gs.value

                if forme == "rect":
                    Asup = w_nbHAs.value * np.pi * (w_phis.value / 1000 / 2) ** 2
                    Ainf = w_nbHAi.value * np.pi * (w_phii.value / 1000 / 2) ** 2
                    section = dict(b=w_b.value, h=w_h.value,
                                   c_inf=w_cinf.value, c_sup=w_csup.value,
                                   As_inf=Ainf, As_sup=Asup,
                                   nb_sup=w_nbHAs.value, nb_inf=w_nbHAi.value)
                else:
                    As_tot = w_nbb.value * np.pi * (w_phib.value / 1000 / 2) ** 2
                    section = dict(D=w_D.value, c_enr=w_enr.value,
                                   nb_barres=w_nbb.value, As_tot=As_tot)

                sollicitations = []
                for ligne in w_sol_txt.value.strip().splitlines():
                    if not ligne.strip():
                        continue
                    parts = [p.strip() for p in ligne.split(",")]
                    if len(parts) >= 2:
                        Ned, Med = float(parts[0]), float(parts[1])
                        lbl = parts[2] if len(parts) >= 3 else f"Cas{len(sollicitations)+1}"
                        sollicitations.append((Ned, Med, lbl))

                (N_arr, M_arr, mat, pts_cles, fibres, arma, Ac, H, zones
                 ) = diagramme_interaction(
                    section_type=forme, section_params=section,
                    fck=fck, fyk=fyk, gamma_c=gc, gamma_s=gs,
                    n_div=400, n_piv_A=200, n_piv_B=300, n_piv_C=150)

                # Mémorise l'état courant pour les appels ultérieurs
                # (ec2.tracer_section(), ec2.tracer(), etc. sans arguments)
                _ETAT_GUI.update(dict(
                    section_type=forme, section_params=section,
                    fck=fck, fyk=fyk, gamma_c=gc, gamma_s=gs,
                    N_arr=N_arr, M_arr=M_arr, mat=mat, pts_cles=pts_cles,
                    fibres=fibres, arma=arma, Ac=Ac, H=H, zones=zones,
                    sollicitations=sollicitations))

                print("═" * 60)
                print("  POINTS CARACTÉRISTIQUES")
                print("═" * 60)
                for k, v in pts_cles.items():
                    print(f"  {k:<42} {v:>10.1f}")
                print("═" * 60)

                fig = tracer_interactif(
                    N_arr, M_arr, mat, pts_cles,
                    section_type=forme, section_params=section,
                    fibres=fibres, arma=arma, Ac=Ac, H=H, zones=zones,
                    sollicitations=sollicitations if sollicitations else None,
                    nom_fichier="diagramme_interactif_gui.html")
                fig.show()

                tracer_section(nom_fichier="section_gui.png")

            except Exception as e:
                print(f"⚠️  Erreur : {type(e).__name__} — {e}")

    btn.on_click(_on_click)

    # ── Boutons d'export (utilisent l'état du dernier "Calculer") ──────────
    btn_pdf  = widgets.Button(description="Export PDF",  icon="file-pdf-o")
    btn_csv  = widgets.Button(description="Export CSV",  icon="table")
    btn_html = widgets.Button(description="Export HTML", icon="globe")

    def _on_export_pdf(_b):
        with out:
            try:
                exporter_pdf(nom_fichier="rapport_gui.pdf")
            except Exception as e:
                print(f"⚠️  Erreur : {type(e).__name__} — {e}")

    def _on_export_csv(_b):
        with out:
            try:
                exporter_csv(nom_fichier="diagramme_gui.csv")
            except Exception as e:
                print(f"⚠️  Erreur : {type(e).__name__} — {e}")

    def _on_export_html(_b):
        with out:
            try:
                fig = tracer_interactif(nom_fichier="diagramme_gui.html")
                fig.show()
            except Exception as e:
                print(f"⚠️  Erreur : {type(e).__name__} — {e}")

    btn_pdf.on_click(_on_export_pdf)
    btn_csv.on_click(_on_export_csv)
    btn_html.on_click(_on_export_html)

    ligne_exports = widgets.HBox([
        widgets.HTML("<b>Exports&nbsp;:</b>&nbsp;"), btn_pdf, btn_csv, btn_html
    ])

    ui = widgets.VBox([
        widgets.HTML("<h3>Diagramme d'interaction N-M — Eurocode 2</h3>"),
        w_forme,
        box_materiaux,
        box_section,
        w_sol_txt,
        btn,
        ligne_exports,
        out,
    ])
    display(ui)


# ═══════════════════════════════════════════════════════════════
# 8.  PROGRAMME PRINCIPAL  —  exemples d'utilisation
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":

    import sys
    # Choisir la section : "rect" ou "circ"  (argument en ligne de commande)
    MODE = sys.argv[1] if len(sys.argv) > 1 else "rect"

    # ══════════════════════════════════════════
    #  SECTION RECTANGULAIRE  (exemple ENPC §13)
    # ══════════════════════════════════════════
    if MODE == "rect":
        fck, fyk = 50.0, 500.0

        phi  = 25e-3           # HA 25
        As   = 3 * np.pi * (phi/2)**2      # 3 HA25 = 14.73 cm²

        sp = dict(
            b      = 0.60,     # m
            h      = 0.60,     # m
            c_inf  = 0.060,    # enrobage axe armature inférieure [m]
            c_sup  = 0.060,    # enrobage axe armature supérieure [m]
            As_inf = As,       # m²
            As_sup = As,       # m²
        )

        sollicitations = [
            ( 3000.0,  800.0, "C1"),    # (N_Ed [kN], M_Ed [kN·m], label)
            ( 2189.0,  842.0, "ENPC"),  # point de l'exemple du cours
            (-400.0,   300.0, "C3"),
        ]

        print("\n  ► MODE : SECTION RECTANGULAIRE (exemple ENPC C50, 60×60, 3HA25/face)")

    # ══════════════════════════════════════════
    #  SECTION CIRCULAIRE
    # ══════════════════════════════════════════
    else:
        fck, fyk = 30.0, 500.0

        phi      = 20e-3           # HA 20
        nb       = 12
        As_tot   = nb * np.pi * (phi/2)**2   # 12 HA20 = 37.70 cm²

        sp = dict(
            D         = 0.60,
            c_enr     = 0.055,
            nb_barres = nb,
            As_tot    = As_tot,
        )

        sollicitations = [
            ( 1500.0, 250.0, "C1"),
            ( 3000.0, 150.0, "C2"),
            ( -300.0, 100.0, "C3"),
        ]

        print("\n  ► MODE : SECTION CIRCULAIRE (Ø60cm, C30, 12HA20)")

    # ── Calcul & tracé ───────────────────────────────────────────────────
    (N_arr, M_arr, mat, pts_cles, fibres, arma, Ac, H, zones) = diagramme_interaction(
        section_type   = MODE,
        section_params = sp,
        fck=fck, fyk=fyk,
        n_div=200, n_piv_A=80, n_piv_B=150, n_piv_C=50,
    )

    tracer(
        N_arr, M_arr, mat, pts_cles,
        section_type   = MODE,
        section_params = sp,
        fibres=fibres, arma=arma, Ac=Ac, H=H,
        zones=zones,
        sollicitations = sollicitations,
        nom_fichier    = f"diagramme_{MODE}_EC2_ENPC.png",
    )

    tracer_interactif(
        N_arr, M_arr, mat, pts_cles,
        section_type   = MODE,
        section_params = sp,
        fibres=fibres, arma=arma, Ac=Ac, H=H,
        zones=zones,
        sollicitations = sollicitations,
        nom_fichier    = f"diagramme_{MODE}_EC2_ENPC_interactif.html",
    )
