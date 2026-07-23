"""
=======================================================================
  JUSTIFICATION À LA FLEXION SIMPLE  —  Eurocode 2 (EN 1992-1-1)
  V1 : section RECTANGULAIRE — ELU, ELS, aciers minimaux
=======================================================================
Réutilise diagramme_interaction_EC2.py comme moteur de capacité ELU
(la flexion simple est le cas particulier N=0 du diagramme d'interaction
N-M déjà construit et validé) — aucune nouvelle formule de résistance
béton/acier n'est réécrite pour l'ELU.

Références normatives (citations exactes extraites de la NF EN 1992-1-1,
texte français, édition consultée) :

  §6.1   Flexion simple et flexion composée — méthode de calcul ELU
         (identique à celle déjà implémentée : parabole-rectangle §3.1.7,
         paliers pivots A/B/C)
  §7.2   Limitation des contraintes ELS :
         σc ≤ k1·fck  (k1 = 0,6, valeur recommandée EC2)  — combinaison caractéristique
         σs ≤ k3·fyk  (k3 = 0,8, valeur recommandée EC2)  — combinaison caractéristique
  §7.3.2 Sections minimales d'armatures (maîtrise de la fissuration) :
         As,min·σs = kc·k·fct,eff·Act                                (7.1)
         kc = 0,4·[1 − rc/(k1·(h/h*)·fct,eff)] ≤ 1  (rc = NEd/(b·h))  (7.2)
  §7.3.4 Calcul de l'ouverture des fissures :
         wk = sr,max·(εsm − εcm)                                     (7.8)
         εsm−εcm = [σs − kt·(fct,eff/ρp,eff)·(1+αe·ρp,eff)] / Es ≥ 0,6·σs/Es   (7.9)
         sr,max = k3·c + k1·k2·k4·φ/ρp,eff                           (7.11)
  §9.2.1.1 Sections minimale et maximale d'armatures (ELU) :
         As,min = 0,26·(fctm/fyk)·bt·d ≥ 0,0013·bt·d                 (9.1N)

  Tableau 3.1 — formules analytiques (et non table interpolée) :
         fcm  = fck + 8
         fctm = 0,30·fck^(2/3)                    si fck ≤ 50 MPa
         fctm = 2,12·ln(1 + fcm/10)                si fck > 50 MPa
         Ecm  = 22 000·(fcm/10)^0,3   [MPa]

  ANNEXE NATIONALE FRANÇAISE (NF EN 1992-1-1/NA, mars 2007) — désormais
  incorporée pour les clauses qui nous concernent :
    §7.2(2)(3)(5)   : k1=0,6 ; k3=0,8 ; k4=1 — confirmés = valeurs EC2 recommandées (inchangé)
    §7.3.1(5)       : Tableau 7.1NF (wmax) — confirmé identique au Tableau 7.1N EC2 (inchangé)
    §7.3.2(2)       : σs = fyk pour le calcul de As,min fissuration — confirmé (inchangé)
    §7.3.4(3) NOTE  : k3 (Expression 7.11) DÉPEND DE L'ENROBAGE au-delà de 25mm
                      (k3 = 3,4 pour c≤25mm, k3 = 3,4·(25/c)^(2/3) au-delà)
                      — CECI DIFFÈRE de la valeur EC2 "recommandée" constante,
                      corrigé dans k3_ANF() ci-dessous.
    §7.3.4(3) NOTE  : la bascule vers l'Expression (7.14) ne se fait que si elle
                      donne un sr,max PLUS GRAND que (7.11) — pas simplement
                      quand l'espacement dépasse 5(c+φ/2) comme le suggère la
                      Figure 7.2 de l'EC2 seul. Corrigé dans ouverture_fissure().
"""

import numpy as np
import diagramme_interaction_EC2 as ec2


# ═══════════════════════════════════════════════════════════════
# 1.  CARACTÉRISTIQUES BÉTON — Tableau 3.1 (formules analytiques)
# ═══════════════════════════════════════════════════════════════

def caracteristiques_beton(fck):
    """fcm, fctm, Ecm selon Tableau 3.1 EC2 (formules exactes, fck en MPa)."""
    fcm = fck + 8.0
    if fck <= 50.0:
        fctm = 0.30 * fck ** (2.0 / 3.0)
    else:
        fctm = 2.12 * np.log(1.0 + fcm / 10.0)
    Ecm = 22000.0 * (fcm / 10.0) ** 0.3   # MPa
    return dict(fcm=fcm, fctm=fctm, Ecm=Ecm)


# ═══════════════════════════════════════════════════════════════
# 2.  VALEURS RECOMMANDÉES EC2 (§7.1N, modifiables si ANF différente)
# ═══════════════════════════════════════════════════════════════

# wmax [mm] — Tableau 7.1N, béton armé, combinaison quasi-permanente
WMAX_TABLE = {
    "X0":  0.4, "XC1": 0.4,
    "XC2": 0.3, "XC3": 0.3, "XC4": 0.3,
    "XD1": 0.3, "XD2": 0.3, "XS1": 0.3, "XS2": 0.3, "XS3": 0.3,
}

COEFFS_RECOMMANDES = dict(
    k1_contrainte_beton=0.6,   # §7.2(2) — limite σc = k1.fck (classes XD/XF/XS)
    k2_fluage=0.45,            # §7.2(3) — linéarité du fluage
    k3_contrainte_acier=0.8,   # §7.2(5) — limite σs = k3.fyk
    kt_court_terme=0.6,        # §7.3.4(2) — chargement courte durée
    kt_long_terme=0.4,         # §7.3.4(2) — chargement longue durée
    k1_adherence=0.8,          # §7.3.4(3) — barres HA
    k2_flexion=0.5,            # §7.3.4(3) — flexion (vs 1,0 en traction pure)
    # k3_espacement : PAS une constante — voir k3_ANF() ci-dessous (ANF §7.3.4(3))
    k4_espacement=0.425,       # §7.3.4(3) — Expression (7.11)
)


def k3_ANF(c_mm):
    """
    Coefficient k3 de l'Expression (7.11), selon l'Annexe Nationale
    française §7.3.4(3) NOTE :
      k3 = 3,4                    si c ≤ 25 mm
      k3 = 3,4 . (25/c)^(2/3)     si c > 25 mm   (c en mm)

    (La valeur EC2 "recommandée" 3,4 constante, quel que soit l'enrobage,
    n'est PAS celle retenue par la France au-delà de 25 mm d'enrobage.)
    """
    if c_mm <= 25.0:
        return 3.4
    return 3.4 * (25.0 / c_mm) ** (2.0 / 3.0)


# ═══════════════════════════════════════════════════════════════
# 2bis.  GÉOMÉTRIE RÉELLE DU FERRAILLAGE (lits multiples)
# ═══════════════════════════════════════════════════════════════

def geometrie_nappe(h, enrobage_nominal, nb_barres_par_lit, phi_mm,
                     entraxe_vertical_mm=None, cote="inf"):
    """
    Calcule la géométrie réelle d'une nappe d'armatures pouvant comporter
    1 ou 2 lits (superposés), à partir des données physiques réelles :
    nombre de barres par lit, diamètre, enrobage nominal (à la génératrice
    de la barre du lit le plus proche du parement).

    entraxe_vertical_mm : distance axe à axe entre lits successifs. Si
    None, valeur par défaut = φ + 20 mm (écartement libre minimal usuel
    entre lits, cf. cours BA / EC2 §8.2 — modifiable si le projet impose
    une valeur différente, p.ex. pour le passage d'un vibreur).

    cote : "inf" (nappe tendue en flexion positive, près du bas) ou
           "sup" (près du haut) — pilote le sens de l'empilement des lits
           (le lit 1 est toujours le plus proche du parement concerné).

    Retourne dict(
        As_total    : aire totale [m²]
        d_eff       : bras de levier réel — distance de la fibre EXTRÊME
                      OPPOSÉE (comprimée) au centroïde de la nappe [m]
                      (c'est le "d" ou "d'" à utiliser dans les calculs)
        c_eff       : enrobage EFFECTIF équivalent au centroïde, mesuré
                      depuis le parement concerné [m] (= h - d_eff si
                      cote="inf", = d_eff... attention, cf. code)
        lits        : liste de dicts {y, nb, phi_mm, As} pour le dessin
                      (y = position réelle en repère centré, + = haut)
    )
    """
    phi = phi_mm / 1000.0
    if entraxe_vertical_mm is None:
        entraxe_mm = phi_mm + 20.0
    else:
        entraxe_mm = entraxe_vertical_mm
    entraxe = entraxe_mm / 1000.0

    lits = []
    for i, nb in enumerate(nb_barres_par_lit):
        if nb <= 0:
            continue
        # distance du CENTRE de la barre au parement concerné
        dist_parement = enrobage_nominal + phi / 2.0 + i * entraxe
        if cote == "inf":
            y = -h / 2.0 + dist_parement
        else:
            y = +h / 2.0 - dist_parement
        As_lit = nb * np.pi * (phi / 2.0) ** 2
        lits.append(dict(y=y, nb=nb, phi_mm=phi_mm, As=As_lit,
                          dist_parement=dist_parement))

    As_total = sum(l["As"] for l in lits)
    if As_total <= 0:
        return dict(As_total=0.0, d_eff=h - enrobage_nominal, c_eff=enrobage_nominal,
                    lits=[])

    y_centroide = sum(l["y"] * l["As"] for l in lits) / As_total

    if cote == "inf":
        d_eff = h / 2.0 - y_centroide          # distance depuis le HAUT (fibre comprimée)
        c_eff = h / 2.0 + y_centroide           # = h - d_eff, depuis le BAS
    else:
        d_eff = h / 2.0 + y_centroide           # distance depuis le BAS
        c_eff = h / 2.0 - y_centroide           # depuis le HAUT

    return dict(As_total=As_total, d_eff=d_eff, c_eff=c_eff, lits=lits)


# ═══════════════════════════════════════════════════════════════
# 3.  CAPACITÉ ELU — réutilise le moteur diagramme_interaction_EC2
# ═══════════════════════════════════════════════════════════════

def capacite_ELU_flexion_simple(section, fck, fyk, gamma_c=1.5, gamma_s=1.15,
                                 n_div=400, n_piv_A=200, n_piv_B=300, n_piv_C=150):
    """
    Moment résistant ELU en flexion simple (N=0), positif et négatif,
    obtenu par interpolation du diagramme d'interaction N-M à N=0.

    section : dict(b, h, c_inf, c_sup, As_inf, As_sup)  [mêmes clés que
              diagramme_interaction_EC2.diagramme_interaction]

    Retourne dict(M_Rd_pos, M_Rd_neg, mat, pts_cles, N_arr, M_arr, zones,
                  fibres, arma, Ac, H)
    """
    (N_arr, M_arr, mat, pts_cles, fibres, arma, Ac, H, zones
     ) = ec2.diagramme_interaction(
        section_type="rect", section_params=section,
        fck=fck, fyk=fyk, gamma_c=gamma_c, gamma_s=gamma_s,
        n_div=n_div, n_piv_A=n_piv_A, n_piv_B=n_piv_B, n_piv_C=n_piv_C)

    # Interpolation des points où N_arr change de signe (contour fermé,
    # généralement 2 croisements : un côté M>0, un côté M<0)
    croisements = []
    for i in range(len(N_arr) - 1):
        if N_arr[i] == 0.0:
            croisements.append(M_arr[i])
        elif N_arr[i] * N_arr[i + 1] < 0:
            t = N_arr[i] / (N_arr[i] - N_arr[i + 1])
            M_interp = M_arr[i] + t * (M_arr[i + 1] - M_arr[i])
            croisements.append(M_interp)

    M_Rd_pos = max([m for m in croisements if m >= 0], default=0.0)
    M_Rd_neg = min([m for m in croisements if m <= 0], default=0.0)

    return dict(M_Rd_pos=M_Rd_pos, M_Rd_neg=M_Rd_neg, mat=mat,
                pts_cles=pts_cles, N_arr=N_arr, M_arr=M_arr, zones=zones,
                fibres=fibres, arma=arma, Ac=Ac, H=H)


def moment_reduit(section, fck, fyk, M_ELU_kNm, gamma_c=1.5, gamma_s=1.15,
                   n_div=400):
    """
    Moment réduit µ = M_Ed/(b.d².fcd) et sa valeur limite µlim (frontière
    pivot B "pur" — acier tendu à εyd pile au moment où le béton atteint
    εcu2, cf. la notion de "pivot B pur" déjà utilisée dans le diagramme
    d'interaction). Si µ > µlim, une section simplement armée (aciers
    tendus seuls) ne suffit plus : des aciers comprimés sont nécessaires.

    µlim est calculé par intégration EXACTE de la loi parabole-rectangle
    (pas par la formule simplifiée du bloc rectangulaire équivalent
    λ=0,8/η=1,0, plus approximative) : on intègre le béton seul (sans
    acier) sur la hauteur comprimée x_lim = d.εcu2/(εcu2+εyd), et on
    ramène le moment de la résultante béton au niveau des aciers tendus.

    Retourne dict(mu, mu_lim, x_lim, d, besoin_aciers_comprimes, M_lim)
    """
    b, h, c_inf = section["b"], section["h"], section["c_inf"]
    d = h - c_inf

    mat = ec2.get_material_params(fck, fyk, gamma_c, gamma_s)
    fcd = mat["fcd"] ;  ecu2 = mat["eps_cu2"] ;  eps_yd = mat["eps_yd"]

    # x_lim : frontière pivot B pur, mesurée depuis la fibre sup. comprimée
    x_lim = d * ecu2 / (ecu2 + eps_yd)

    # Intégration exacte (béton seul, sans acier) sur H = h, en extrapolant
    # linéairement le champ de déformation pour que ε=0 exactement à y=x_lim
    # sous la fibre sup. (convention H centrée du moteur diagramme_interaction)
    fibres = ec2.fibres_rect(b, h, n_div)
    eps_top = ecu2
    pente = ecu2 / x_lim
    eps_bot = ecu2 - pente * h
    N_lim, M_lim_centre = ec2.compute_NM(eps_top, eps_bot, h, fibres, [], mat)

    # Report du moment du centre de section vers le niveau des aciers tendus
    # (Varignon : M/A = M/G - N.(y_A - y_G) ; ici y_A = h/2 - d, position des
    # aciers tendus sous le centre, DONC NÉGATIVE dans la convention centrée)
    y_As = h / 2.0 - d
    M_lim = M_lim_centre - N_lim * y_As   # kN.m, moment résistant béton seul / aciers tendus

    mu_lim = (M_lim * 1e3) / (b * d ** 2 * fcd * 1e6)
    mu = (abs(M_ELU_kNm) * 1e3) / (b * d ** 2 * fcd * 1e6)

    return dict(mu=mu, mu_lim=mu_lim, x_lim=x_lim, d=d, M_lim=M_lim,
                besoin_aciers_comprimes=mu > mu_lim)


def _moment_beton_seul(x, h, fibres, mat, d):
    """
    Moment résistant du béton seul (sans acier), ramené au niveau des
    aciers tendus, pour une profondeur d'axe neutre x donnée (pivot B :
    εc,sup=εcu2 fixé). Fonction interne à moment_reduit()/acier_theorique_ELU().
    Retourne (Fc [kN], M_A [kN.m], eps_bot_extreme).
    """
    ecu2 = mat["eps_cu2"]
    pente = ecu2 / x
    eps_bot = ecu2 - pente * h
    N_c, M_c_centre = ec2.compute_NM(ecu2, eps_bot, h, fibres, [], mat)
    y_As = h / 2.0 - d
    M_A = M_c_centre - N_c * y_As
    return N_c, M_A, eps_bot


def acier_theorique_ELU(section, fck, fyk, M_ELU_kNm, gamma_c=1.5, gamma_s=1.15,
                         n_div=400, tol=1e-6, max_iter=60):
    """
    Calcule la section d'acier théoriquement NÉCESSAIRE pour équilibrer
    M_ELU_kNm (dimensionnement, par opposition à la vérification d'un
    ferraillage déjà choisi) : résolution exacte (parabole-rectangle,
    pas le bloc simplifié) par recherche de la profondeur d'axe neutre x
    telle que le moment résistant du béton seul, ramené au niveau des
    aciers tendus, égale M_Ed.

    Si µ > µlim (aciers comprimés nécessaires, cf. moment_reduit), la
    section est calculée en double armature selon la méthode classique :
    - As1 équilibre M_lim (béton seul à x_lim) avec l'acier tendu à fyd
    - le complément ΔM = M_Ed - M_lim est équilibré par le couple
      (aciers comprimés As', aciers tendus additionnels As2) au bras de
      levier (d-d')

    Retourne dict(As_theorique, As1, As2, As_comprime_theorique, x, cas)
    où cas ∈ {"simple", "double"}.
    """
    b, h, c_inf = section["b"], section["h"], section["c_inf"]
    c_sup = section["c_sup"]
    d = h - c_inf
    d_prime = c_sup

    mat = ec2.get_material_params(fck, fyk, gamma_c, gamma_s)
    fyd = mat["fyd"] ;  ecu2 = mat["eps_cu2"] ;  eps_yd = mat["eps_yd"]
    fibres = ec2.fibres_rect(b, h, n_div)
    M_Ed = abs(M_ELU_kNm)

    x_lim = d * ecu2 / (ecu2 + eps_yd)
    _, M_lim, _ = _moment_beton_seul(x_lim, h, fibres, mat, d)

    if M_Ed <= M_lim:
        # ── Cas simple armature : recherche de x par dichotomie ────────
        lo, hi = 1e-4, x_lim
        _, M_lo, _ = _moment_beton_seul(lo, h, fibres, mat, d)
        _, M_hi, _ = _moment_beton_seul(hi, h, fibres, mat, d)
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            _, M_mid, _ = _moment_beton_seul(mid, h, fibres, mat, d)
            if abs(M_mid - M_Ed) < tol * max(1.0, M_Ed):
                break
            if M_mid > M_Ed:
                hi = mid
            else:
                lo = mid
        x = mid
        Fc, M_x, eps_bot_extreme = _moment_beton_seul(x, h, fibres, mat, d)
        # déformation de l'acier tendu (niveau y = -h/2+c_inf)
        pente = ecu2 / x
        eps_s = ecu2 - pente * (h - c_inf)   # = eps au niveau des aciers (>0 conv. compression)
        sigma_s_t = abs(ec2.sigma_s(eps_s, mat))
        As_theo = Fc / (sigma_s_t * 1e3)   # Fc[kN] / (sigma[MPa]=1e3 kN/m²) -> m²

        return dict(As_theorique=As_theo, As1=As_theo, As2=0.0,
                    As_comprime_theorique=0.0, x=x, d=d, cas="simple")

    else:
        # ── Cas double armature ─────────────────────────────────────────
        Fc_lim, _, eps_bot_lim = _moment_beton_seul(x_lim, h, fibres, mat, d)
        As1 = Fc_lim / (fyd * 1e3)   # tendu à fyd (pivot B pur -> juste à εyd)

        # déformation des aciers comprimés à x=x_lim (profondeur d' sous la
        # fibre sup. comprimée) : eps(y_from_top) = ecu2.(1 - y_from_top/x_lim)
        pente = ecu2 / x_lim
        eps_sc = ecu2 - pente * d_prime
        sigma_sc = ec2.sigma_s(eps_sc, mat)  # compression -> positif

        delta_M = M_Ed - M_lim
        levier = d - d_prime
        As_comprime = delta_M / (sigma_sc * 1e3 * levier) if levier > 0 else np.inf
        As2 = delta_M / (fyd * 1e3 * levier) if levier > 0 else np.inf
        As_theo = As1 + As2

        return dict(As_theorique=As_theo, As1=As1, As2=As2,
                    As_comprime_theorique=As_comprime, x=x_lim, d=d, cas="double")


def etat_deformation_ELU(section, fck, fyk, gamma_c=1.5, gamma_s=1.15,
                          M_ELU_kNm=None, n_div=400, tol=1e-6, max_iter=60):
    """
    Détermine l'état de déformation réel à l'ELU pour le ferraillage donné
    (vérification), en supposant le pivot B actif (fibre extrême comprimée
    à εcu2 — hypothèse standard pour une section normalement armée ; si le
    solveur ne trouve pas de racine dans ce domaine, la section est très
    probablement pilotée par le pivot A — cf. valeur de retour None).

    M_ELU_kNm sert uniquement à choisir le signe (quelle nappe est tendue) ;
    l'état renvoyé est l'état À LA RUINE (M=M_Rd côté correspondant), pas
    l'état sous M_Ed — c'est l'état conventionnellement montré sur un
    "diagramme des déformations" (cf. cahier des charges).

    Retourne dict(eps_top, eps_bot, x, eps_sc, eps_st, y_sup, y_inf, N, M)
    ou None si aucune racine trouvée dans le domaine pivot B.
    """
    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]
    As_inf, As_sup = section["As_inf"], section["As_sup"]

    mat = ec2.get_material_params(fck, fyk, gamma_c, gamma_s)
    fibres = ec2.fibres_rect(b, h, n_div)
    arma = ec2.armatures_rect(h, c_inf, c_sup, As_inf, As_sup)
    ecu2 = mat["eps_cu2"] ;  eud = mat["eps_ud"]

    positif = (M_ELU_kNm is None) or (M_ELU_kNm >= 0)

    def N_de(eps_var):
        if positif:
            N, _ = ec2.compute_NM(ecu2, eps_var, h, fibres, arma, mat)
        else:
            N, _ = ec2.compute_NM(eps_var, ecu2, h, fibres, arma, mat)
        return N

    lo, hi = -eud, ecu2
    N_lo, N_hi = N_de(lo), N_de(hi)
    if N_lo > 0 or N_hi < 0:
        return None  # pas de racine en domaine pivot B -> pivot A probable

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if abs(N_de(mid)) < tol * max(1.0, abs(N_lo)):
            break
        if N_de(mid) > 0:
            hi = mid
        else:
            lo = mid

    eps_top, eps_bot = (ecu2, mid) if positif else (mid, ecu2)
    N_eq, M_eq = ec2.compute_NM(eps_top, eps_bot, h, fibres, arma, mat)

    x = ecu2 * h / (ecu2 - eps_bot) if positif else ecu2 * h / (ecu2 - eps_top)

    def eps_a(y):
        return eps_top + (eps_bot - eps_top) * (h / 2 - y) / h

    y_sup = h / 2 - c_sup
    y_inf = -h / 2 + c_inf

    return dict(eps_top=eps_top, eps_bot=eps_bot, x=x, N=N_eq, M=M_eq,
                y_sup=y_sup, y_inf=y_inf, eps_sup=eps_a(y_sup),
                eps_inf=eps_a(y_inf), h=h, positif=positif)


# ═══════════════════════════════════════════════════════════════
# 4.  ANALYSE ÉLASTIQUE SECTION FISSURÉE — pour l'ELS (§7.2, §7.3.4)
# ═══════════════════════════════════════════════════════════════

def axe_neutre_elastique_fissure(b, d, d_prime, As, As_prime, n):
    """
    Position x de l'axe neutre élastique (section fissurée, béton tendu
    négligé), méthode des sections homogénéisées classique :

        b·x²/2 + (n-1)·As'·(x-d') − n·As·(d-x) = 0

    d      : hauteur utile des aciers tendus (= h - c_inf)
    d_prime: distance fibre comprimée -> aciers comprimés (= c_sup)
    As, As_prime : aciers tendus, aciers comprimés [m²]
    n      : coefficient d'équivalence Es/Ecm

    Résolution par équation du 2nd degré : A·x² + B·x + C = 0
    """
    A = b / 2.0
    B = (n - 1.0) * As_prime + n * As
    C = -((n - 1.0) * As_prime * d_prime + n * As * d)
    disc = B ** 2 - 4 * A * C
    x = (-B + np.sqrt(disc)) / (2 * A)
    return x


def inertie_fissuree(b, x, d, d_prime, As, As_prime, n):
    """Moment d'inertie de la section homogénéisée fissurée / axe neutre."""
    I_beton = b * x ** 3 / 3.0
    I_ast   = n * As * (d - x) ** 2
    I_asc   = (n - 1.0) * As_prime * (x - d_prime) ** 2 if As_prime > 0 else 0.0
    return I_beton + I_ast + I_asc


def contraintes_ELS(section, fck, M_ELS_kNm, Es=200_000.0, n_impose=None):
    """
    Contraintes élastiques en section fissurée sous moment de service
    M_ELS (combinaison caractéristique, en kN·m — signe : positif = fibre
    sup. comprimée, aciers inf. tendus, comme pour le diagramme d'inter-
    action).

    n_impose : coefficient d'équivalence n=Es/Ec à utiliser directement
    (p.ex. 15, valeur forfaitaire usuelle — cf. cahier des charges). Si
    None (défaut), n est calculé précisément à partir de Ecm(fck) selon
    le Tableau 3.1 EC2 (méthode plus rigoureuse mais moins "standard"
    que le n=15 traditionnellement utilisé en préconception).

    Retourne dict(x, I_fiss, sigma_c, sigma_st, sigma_sc, n, Ecm)
    """
    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]
    As_inf, As_sup = section["As_inf"], section["As_sup"]
    d = h - c_inf

    carac = caracteristiques_beton(fck)
    Ecm = carac["Ecm"]
    n = n_impose if n_impose is not None else Es / Ecm

    M = abs(M_ELS_kNm) * 1e3  # kN.m -> N.m ... on travaille en unités SI (N, m)
    # NB : b, h en m -> x en m ; on garde M en N.m, b en m => contraintes en Pa,
    # reconverties en MPa à la fin.

    if M_ELS_kNm >= 0:
        As_t, As_c, c_t, c_c = As_inf, As_sup, c_inf, c_sup
    else:
        As_t, As_c, c_t, c_c = As_sup, As_inf, c_sup, c_inf

    x = axe_neutre_elastique_fissure(b, d, c_c, As_t, As_c, n)
    I = inertie_fissuree(b, x, d, c_c, As_t, As_c, n)

    sigma_c  = (M * x) / I / 1e6          # MPa, fibre la plus comprimée
    sigma_st = n * M * (d - x) / I / 1e6  # MPa, acier tendu
    sigma_sc = ((n - 1.0) * M * (x - c_c) / I / 1e6) if As_c > 0 else 0.0

    return dict(x=x, I_fiss=I, sigma_c=sigma_c, sigma_st=sigma_st,
                sigma_sc=sigma_sc, n=n, Ecm=Ecm, d=d)


def acier_theorique_ELS(section, fck, fyk, M_ELS_kNm, gamma_c=1.5, gamma_s=1.15,
                         Es=200_000.0, n_impose=None, tol=1e-5, max_iter=60):
    """
    Calcule la section d'acier tendu théoriquement NÉCESSAIRE pour que
    les DEUX limites de contrainte ELS (§7.2) soient respectées :
      σc ≤ k1.fck   et   σs ≤ k3.fyk

    Méthode : pour chacun des deux critères pris séparément, on cherche
    par dichotomie l'aire d'acier tendu As qui amène exactement la
    contrainte concernée à sa valeur limite (σs(As) est décroissante en
    As, σc(As) croît légèrement puis se stabilise) ; l'As théorique
    retenu est le plus grand des deux (celui qui gouverne), garantissant
    que les deux critères sont satisfaits simultanément.

    Retourne dict(As_theorique_ELS, As_critere_acier, As_critere_beton,
                  critere_gouvernant, sigma_c_lim, sigma_s_lim)
    """
    positif = M_ELS_kNm >= 0
    sigma_c_lim = COEFFS_RECOMMANDES["k1_contrainte_beton"] * fck
    sigma_s_lim = COEFFS_RECOMMANDES["k3_contrainte_acier"] * fyk

    def _section_avec_As(As_trial):
        sec2 = dict(section)
        if positif:
            sec2["As_inf"] = As_trial
        else:
            sec2["As_sup"] = As_trial
        return sec2

    def _contraintes(As_trial):
        r = contraintes_ELS(_section_avec_As(As_trial), fck, M_ELS_kNm,
                            Es=Es, n_impose=n_impose)
        return r["sigma_c"], r["sigma_st"]

    def _bissection(critere_fn, lo=1e-5, hi=0.05):
        # critere_fn(As) doit être décroissante et changer de signe sur [lo,hi]
        f_lo, f_hi = critere_fn(lo), critere_fn(hi)
        if f_lo < 0:
            # même l'aire minimale suffit -> pas de contrainte active, on
            # renvoie une valeur très faible (section non gouvernante)
            return lo
        if f_hi > 0:
            # même une aire énorme ne suffit pas (cas très rare/dégénéré)
            hi = 0.5  # dernier recours, section manifestement énorme
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            f_mid = critere_fn(mid)
            if abs(f_mid) < tol:
                break
            if f_mid > 0:
                lo = mid
            else:
                hi = mid
        return mid

    As_acier = _bissection(lambda As: _contraintes(As)[1] - sigma_s_lim)
    As_beton = _bissection(lambda As: _contraintes(As)[0] - sigma_c_lim)

    As_theorique_ELS = max(As_acier, As_beton)
    critere_gouvernant = "acier (σs)" if As_acier >= As_beton else "béton (σc)"

    return dict(As_theorique_ELS=As_theorique_ELS, As_critere_acier=As_acier,
                As_critere_beton=As_beton, critere_gouvernant=critere_gouvernant,
                sigma_c_lim=sigma_c_lim, sigma_s_lim=sigma_s_lim)


# ═══════════════════════════════════════════════════════════════
# 5.  ACIERS MINIMAUX
# ═══════════════════════════════════════════════════════════════

def acier_min_ELU(section, fck, fyk):
    """§9.2.1.1 (9.1N) : As,min = max(0,26.fctm/fyk.bt.d ; 0,0013.bt.d)."""
    b, h, c_inf = section["b"], section["h"], section["c_inf"]
    d = h - c_inf
    fctm = caracteristiques_beton(fck)["fctm"]
    As_min_1 = 0.26 * (fctm / fyk) * b * d
    As_min_2 = 0.0013 * b * d
    return dict(As_min=max(As_min_1, As_min_2), As_min_1=As_min_1, As_min_2=As_min_2, d=d)


def acier_min_fissuration(section, fck, fyk, k=1.0, k1_contrainte_beton=0.6):
    """
    §7.3.2 (7.1)-(7.2), cas flexion simple (NEd=0 -> rc=0 -> kc=0,4),
    section rectangulaire, Act pris égal à b.h/2 (zone tendue de la
    section brute non fissurée en flexion simple symétrique élastique).

    rs pris égal à fyk (valeur recommandée en l'absence de limitation
    d'ouverture de fissure imposant une contrainte inférieure — voir
    ouverture_fissure() pour la vérification complète wk).
    """
    b, h = section["b"], section["h"]
    fctm = caracteristiques_beton(fck)["fctm"]
    Act = b * h / 2.0
    kc = 0.4  # NEd=0 -> rc=0 (flexion simple, cf §7.3.2 (7.2))
    rs = fyk
    As_min = kc * k * fctm * Act / rs
    return dict(As_min=As_min, Act=Act, kc=kc, k=k, fctm=fctm, rs=rs)


# ═══════════════════════════════════════════════════════════════
# 6.  OUVERTURE DE FISSURE — §7.3.4
# ═══════════════════════════════════════════════════════════════

def ouverture_fissure(section, fck, M_ELS_kNm, nb_barres_tendues, phi_barre_mm,
                       duree_chargement="long_terme", Es=200_000.0,
                       coeffs=None, n_impose=None):
    """
    Calcule wk selon (7.8)-(7.14).

    nb_barres_tendues, phi_barre_mm : nombre et diamètre des barres de la
    nappe tendue sous M_ELS (utilisés pour φ et l'espacement).
    duree_chargement : "court_terme" (kt=0,6) ou "long_terme" (kt=0,4).
    n_impose : coefficient d'équivalence n=Es/Ec forfaitaire (cf.
    contraintes_ELS) ; None = calcul précis via Ecm(fck).
    """
    if coeffs is None:
        coeffs = COEFFS_RECOMMANDES

    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]

    carac = caracteristiques_beton(fck)
    fct_eff = carac["fctm"]

    els = contraintes_ELS(section, fck, M_ELS_kNm, Es=Es, n_impose=n_impose)
    x = els["x"] ;  d = els["d"] ;  n_eq = els["n"]
    sigma_s = els["sigma_st"]

    c_tendu = c_inf if M_ELS_kNm >= 0 else c_sup
    phi = phi_barre_mm / 1000.0  # m

    # Ac,eff : hc,ef = min(2,5(h-d), (h-x)/3, h/2)
    hc_ef = min(2.5 * (h - d), (h - x) / 3.0, h / 2.0)
    Ac_eff = b * hc_ef

    As_t = nb_barres_tendues * np.pi * (phi / 2) ** 2
    rho_p_eff = As_t / Ac_eff if Ac_eff > 0 else np.inf

    kt = coeffs["kt_court_terme"] if duree_chargement == "court_terme" else coeffs["kt_long_terme"]
    alpha_e = n_eq

    esm_ecm = (sigma_s - kt * (fct_eff / rho_p_eff) * (1 + alpha_e * rho_p_eff)) / Es
    esm_ecm = max(esm_ecm, 0.6 * sigma_s / Es)   # borne basse (7.9)

    # Espacement des barres (entraxe) — indicatif, conservé pour le rapport,
    # mais NE PILOTE PLUS le choix de formule (cf. correction ANF ci-dessous)
    if nb_barres_tendues > 1:
        entraxe = (b - 2 * c_tendu - phi) / (nb_barres_tendues - 1)
    else:
        entraxe = b - 2 * c_tendu
    limite_entraxe = 5 * (c_tendu + phi / 2)

    # sr,max — Expression (7.11), avec k3 dépendant de l'enrobage (ANF §7.3.4(3))
    k3 = k3_ANF(c_tendu * 1000)  # c en mm pour la formule ANF
    sr_max_7_11 = (k3 * c_tendu
                   + coeffs["k1_adherence"] * coeffs["k2_flexion"]
                   * coeffs["k4_espacement"] * phi / rho_p_eff)

    # sr,max — Expression (7.14), enveloppe supérieure
    sr_max_7_14 = 1.3 * (h - x)

    # Règle ANF §7.3.4(3) NOTE : l'Expression (7.14) ne s'applique QUE si elle
    # donne une valeur PLUS GRANDE que (7.11) — sinon (7.11) reste applicable
    # même si l'espacement dépasse 5(c+φ/2). On ne bascule donc plus sur un
    # critère d'espacement, mais sur la comparaison directe des deux valeurs.
    if sr_max_7_14 > sr_max_7_11:
        sr_max = sr_max_7_14
        formule = "7.14 (enveloppe, > 7.11 — ANF §7.3.4(3))"
    else:
        sr_max = sr_max_7_11
        formule = f"7.11 (k3={k3:.3f}{' [ANF, c>25mm]' if c_tendu*1000 > 25 else ''})"

    wk = sr_max * esm_ecm

    return dict(wk=wk * 1000, sr_max=sr_max * 1000, esm_ecm=esm_ecm,
                sigma_s=sigma_s, rho_p_eff=rho_p_eff, Ac_eff=Ac_eff,
                hc_ef=hc_ef, x=x, formule=formule, entraxe=entraxe * 1000,
                limite_entraxe=limite_entraxe * 1000,
                sr_max_7_11=sr_max_7_11 * 1000, sr_max_7_14=sr_max_7_14 * 1000,
                k3=k3)


# ═══════════════════════════════════════════════════════════════
# 7.  ORCHESTRATEUR — justification complète
# ═══════════════════════════════════════════════════════════════

def justifier_flexion_simple(section, fck, fyk,
                              M_ELU_kNm, M_ELS_kNm,
                              nb_barres_tendues, phi_barre_mm,
                              classe_exposition="XC1",
                              gamma_c=1.5, gamma_s=1.15,
                              duree_chargement="long_terme",
                              n_impose=None):
    """
    Justification complète d'une section rectangulaire en flexion simple :
    moment réduit µ / besoin d'aciers comprimés, capacité ELU, état de
    déformation ELU (diagramme), contraintes ELS, aciers théoriques et
    minimaux (ELU + ELS/fissuration), ouverture de fissure wk vs
    wmax(classe d'exposition).

    n_impose : coefficient d'équivalence acier/béton n=Es/Ec à utiliser
    pour les calculs ELS (contraintes, ouverture de fissure, As théorique
    ELS). None (défaut) = calcul précis n=Es/Ecm(fck) ; une valeur
    forfaitaire usuelle est n=15 (cf. cahier des charges).

    Retourne un dict structuré avec tous les résultats + verdicts booléens.
    """
    resultats = {}

    # --- Moment réduit / besoin d'aciers comprimés ---
    mr = moment_reduit(section, fck, fyk, M_ELU_kNm, gamma_c, gamma_s)
    resultats["moment_reduit"] = mr

    # --- Acier théorique nécessaire à l'ELU (dimensionnement) ---
    resultats["acier_theorique"] = acier_theorique_ELU(section, fck, fyk, M_ELU_kNm, gamma_c, gamma_s)

    # --- Acier théorique nécessaire à l'ELS (dimensionnement, §7.2) ---
    resultats["acier_theorique_ELS"] = acier_theorique_ELS(
        section, fck, fyk, M_ELS_kNm, gamma_c, gamma_s, n_impose=n_impose)

    # --- État de déformation ELU (pour le diagramme) ---
    etat = etat_deformation_ELU(section, fck, fyk, gamma_c, gamma_s, M_ELU_kNm)
    resultats["deformation_ELU"] = etat  # peut être None si pivot A gouverne

    # --- ELU (capacité) ---
    elu = capacite_ELU_flexion_simple(section, fck, fyk, gamma_c, gamma_s)
    M_Rd = elu["M_Rd_pos"] if M_ELU_kNm >= 0 else elu["M_Rd_neg"]
    resultats["ELU"] = dict(
        M_Ed=M_ELU_kNm, M_Rd=M_Rd,
        verifie=abs(M_ELU_kNm) <= abs(M_Rd),
        taux=abs(M_ELU_kNm) / abs(M_Rd) if M_Rd != 0 else np.inf,
        detail=elu)

    # --- ELS : contraintes ---
    els = contraintes_ELS(section, fck, M_ELS_kNm, n_impose=n_impose)
    k1c = COEFFS_RECOMMANDES["k1_contrainte_beton"]
    k3s = COEFFS_RECOMMANDES["k3_contrainte_acier"]
    sigma_c_lim = k1c * fck
    sigma_s_lim = k3s * fyk
    resultats["ELS_contraintes"] = dict(
        sigma_c=els["sigma_c"], sigma_c_lim=sigma_c_lim,
        verifie_beton=els["sigma_c"] <= sigma_c_lim,
        sigma_s=els["sigma_st"], sigma_s_lim=sigma_s_lim,
        verifie_acier=els["sigma_st"] <= sigma_s_lim,
        detail=els)

    # --- Aciers minimaux (ELU §9.2.1.1, et ELS/fissuration §7.3.2) ---
    amin_elu = acier_min_ELU(section, fck, fyk)
    amin_els = acier_min_fissuration(section, fck, fyk)
    As_inf = section["As_inf"] ;  As_sup = section["As_sup"]
    As_tendu = As_inf if M_ELU_kNm >= 0 else As_sup
    resultats["aciers_minimaux"] = dict(
        As_tendu_reel=As_tendu,
        As_min_ELU=amin_elu["As_min"], verifie_ELU=As_tendu >= amin_elu["As_min"],
        As_min_ELS=amin_els["As_min"], verifie_ELS=As_tendu >= amin_els["As_min"],
        # alias conservés pour compatibilité avec du code appelant existant :
        As_min_fissuration=amin_els["As_min"], verifie_fissuration=As_tendu >= amin_els["As_min"],
        detail_ELU=amin_elu, detail_ELS=amin_els)

    # --- Ouverture de fissure ---
    wmax = WMAX_TABLE.get(classe_exposition, 0.3)
    fiss = ouverture_fissure(section, fck, M_ELS_kNm, nb_barres_tendues,
                              phi_barre_mm, duree_chargement=duree_chargement,
                              n_impose=n_impose)
    resultats["fissuration"] = dict(
        wk=fiss["wk"], wmax=wmax, verifie=fiss["wk"] <= wmax,
        classe_exposition=classe_exposition, detail=fiss)

    resultats["verifie_global"] = all([
        resultats["ELU"]["verifie"],
        resultats["ELS_contraintes"]["verifie_beton"],
        resultats["ELS_contraintes"]["verifie_acier"],
        resultats["aciers_minimaux"]["verifie_ELU"],
        resultats["aciers_minimaux"]["verifie_ELS"],
        resultats["fissuration"]["verifie"],
        not resultats["moment_reduit"]["besoin_aciers_comprimes"],
    ])

    return resultats


def imprimer_rapport(resultats):
    """Affiche un rapport texte lisible du dict retourné par justifier_flexion_simple()."""
    def tag(ok):
        return "✔ VÉRIFIÉ    " if ok else "✘ NON VÉRIFIÉ"

    print("═" * 68)
    print("  JUSTIFICATION FLEXION SIMPLE — EC2")
    print("═" * 68)

    mr = resultats["moment_reduit"]
    print(f"\n[Moment réduit]  µ = {mr['mu']:.4f}   µlim = {mr['mu_lim']:.4f}"
          f"   {'✘ ACIERS COMPRIMÉS NÉCESSAIRES' if mr['besoin_aciers_comprimes'] else '✔ Section simplement armée suffisante'}")

    at = resultats["acier_theorique"]
    print(f"\n[Acier théorique nécessaire à l'ELU]  (dimensionnement, cas {at['cas']})")
    print(f"  As théorique tendu = {at['As_theorique']*1e4:.2f} cm²")
    if at["cas"] == "double":
        print(f"    dont As1 (équilibre M_lim) = {at['As1']*1e4:.2f} cm², "
              f"As2 (complément) = {at['As2']*1e4:.2f} cm²")
        print(f"  As théorique comprimé       = {at['As_comprime_theorique']*1e4:.2f} cm²")

    etat = resultats["deformation_ELU"]
    if etat is not None:
        print(f"\n[État de déformation à l'ELU]  (pivot B, à la ruine)")
        print(f"  x (axe neutre / fibre comprimée) = {etat['x']*1000:.1f} mm")
        print(f"  εsup = {etat['eps_sup']*1e3:+.2f}‰   εinf = {etat['eps_inf']*1e3:+.2f}‰")
    else:
        print(f"\n[État de déformation à l'ELU]  pivot A probable (section peu armée) — non calculé")

    e = resultats["ELU"]
    print(f"\n[ELU]  M_Ed = {e['M_Ed']:8.1f} kN·m   M_Rd = {e['M_Rd']:8.1f} kN·m"
          f"   (taux {e['taux']*100:5.1f}%)   {tag(e['verifie'])}")

    s = resultats["ELS_contraintes"]
    print(f"\n[ELS - contraintes]")
    print(f"  σc = {s['sigma_c']:6.2f} MPa  ≤  {s['sigma_c_lim']:5.2f} MPa (k1.fck)   {tag(s['verifie_beton'])}")
    print(f"  σs = {s['sigma_s']:6.1f} MPa  ≤  {s['sigma_s_lim']:5.1f} MPa (k3.fyk)   {tag(s['verifie_acier'])}")

    a = resultats["aciers_minimaux"]
    print(f"\n[Aciers minimaux]  As (nappe tendue) = {a['As_tendu_reel']*1e4:.2f} cm²")
    print(f"  As,min ELU (§9.2.1.1)          = {a['As_min_ELU']*1e4:6.2f} cm²   {tag(a['verifie_ELU'])}")
    print(f"  As,min ELS/fissuration (§7.3.2) = {a['As_min_ELS']*1e4:6.2f} cm²   {tag(a['verifie_ELS'])}")

    f = resultats["fissuration"]
    print(f"\n[Ouverture de fissure]  classe {f['classe_exposition']}")
    print(f"  wk = {f['wk']:.3f} mm  ≤  wmax = {f['wmax']:.2f} mm   {tag(f['verifie'])}")
    print(f"  (formule : {f['detail']['formule']})")

    print("\n" + "═" * 68)
    verdict = "✔  SECTION JUSTIFIÉE" if resultats["verifie_global"] else "✘  SECTION NON JUSTIFIÉE"
    print(f"  {verdict}")
    print("═" * 68)


# ═══════════════════════════════════════════════════════════════
# 7bis.  SCHÉMA DE SECTION DÉTAILLÉ — lits réels, enrobage, bras de levier
# ═══════════════════════════════════════════════════════════════

def schema_section_detaille(b, h, geom_inf=None, geom_sup=None, nom_fichier=None):
    """
    Schéma de section coté, à l'échelle, montrant le ferraillage RÉEL :
    - chaque lit dessiné séparément (bon nombre de barres, bon diamètre)
    - enrobage réel coté (jusqu'au nu de la barre du 1er lit)
    - bras de levier réel d_eff coté (jusqu'au centroïde de la nappe,
      qui ne coïncide avec le 1er lit que s'il n'y a qu'un seul lit)

    geom_inf, geom_sup : dict renvoyés par geometrie_nappe() (ou None si
    la nappe est vide).

    Retourne la figure matplotlib (fig).
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    sc = 1.0 / h
    bS, hS = b * sc, h * sc

    fig, ax = plt.subplots(figsize=(5.5, 8.5))
    ax.add_patch(mpatches.Rectangle((-bS/2, 0), bS, hS, fc="#eef1f7", ec="#333", lw=1.6))

    r_scale = min(bS, hS)

    def dessiner_nappe(geom, cote):
        if geom is None or geom["As_total"] <= 0:
            return
        for i, lit in enumerate(geom["lits"]):
            y_S = (lit["y"] + h/2) * sc   # repère section (0 en bas)
            r_b = 0.018 + 0.0035 * lit["phi_mm"] / 12.0
            nb = lit["nb"]
            marge = bS * 0.08
            xs = np.linspace(-bS/2 + marge, bS/2 - marge, nb) if nb > 1 else [0.0]
            for xb in xs:
                ax.add_patch(plt.Circle((xb, y_S), r_b, color="#b71c1c", zorder=5))
            label = f"Lit {i+1} : {nb}HA{int(lit['phi_mm'])}"
            side = -1 if cote == "inf" else 1
            ax.annotate(label, (bS/2 + marge, y_S),
                        xytext=(bS/2 + 0.10, y_S), fontsize=8.2,
                        va="center", ha="left", color="#7a1414")

        # Enrobage réel (jusqu'au nu du lit 1)
        lit1 = geom["lits"][0]
        y1_S = (lit1["y"] + h/2) * sc
        phi1 = lit1["phi_mm"] / 1000.0 * sc
        if cote == "inf":
            y_nu = 0.0
            y_bar_nu = y1_S - phi1/2
            ax.annotate("", xy=(-bS/2 - 0.14, y_nu), xytext=(-bS/2 - 0.14, y_bar_nu),
                        arrowprops=dict(arrowstyle="<->", color="#2E7D32", lw=1))
            ax.text(-bS/2 - 0.17, (y_nu+y_bar_nu)/2, f"c={geom['lits'][0]['dist_parement']*1000 - lit1['phi_mm']/2:.0f}mm",
                    fontsize=7.8, color="#2E7D32", ha="right", va="center")
        else:
            y_nu = hS
            y_bar_nu = y1_S + phi1/2
            ax.annotate("", xy=(-bS/2 - 0.14, y_nu), xytext=(-bS/2 - 0.14, y_bar_nu),
                        arrowprops=dict(arrowstyle="<->", color="#2E7D32", lw=1))
            ax.text(-bS/2 - 0.17, (y_nu+y_bar_nu)/2, f"c={geom['lits'][0]['dist_parement']*1000 - lit1['phi_mm']/2:.0f}mm",
                    fontsize=7.8, color="#2E7D32", ha="right", va="center")

        # Bras de levier réel (jusqu'au centroïde de la nappe)
        if cote == "inf":
            y_centroid_S = (h - geom["d_eff"]) * sc   # centroïde, repère bas=0
        else:
            y_centroid_S = geom["d_eff"] * sc
        if cote == "inf":
            y_opp = hS
            ax.annotate("", xy=(bS/2 + 0.42, y_opp), xytext=(bS/2 + 0.42, y_centroid_S),
                        arrowprops=dict(arrowstyle="<->", color="#1565C0", lw=1.1))
            ax.text(bS/2 + 0.45, (y_opp+y_centroid_S)/2, f"d={geom['d_eff']*1000:.0f}mm",
                    fontsize=8.5, color="#1565C0", va="center", fontweight="bold")
        else:
            y_opp = 0.0
            ax.annotate("", xy=(bS/2 + 0.42, y_opp), xytext=(bS/2 + 0.42, y_centroid_S),
                        arrowprops=dict(arrowstyle="<->", color="#1565C0", lw=1.1))
            ax.text(bS/2 + 0.45, (y_opp+y_centroid_S)/2, f"d'={geom['d_eff']*1000:.0f}mm",
                    fontsize=8.5, color="#1565C0", va="center", fontweight="bold")
        # marqueur du centroïde (si >1 lit, distinct des barres)
        if len(geom["lits"]) > 1:
            ax.plot(0, y_centroid_S, "+", color="#1565C0", ms=12, mew=2, zorder=6)

    dessiner_nappe(geom_inf, "inf")
    dessiner_nappe(geom_sup, "sup")

    # Cotes b, h
    ax.annotate("", xy=(-bS/2, -0.08), xytext=(bS/2, -0.08),
                arrowprops=dict(arrowstyle="<->", color="#333", lw=1))
    ax.text(0, -0.13, f"b={b*100:.0f}cm", ha="center", fontsize=9)
    ax.annotate("", xy=(bS/2+0.75, 0), xytext=(bS/2+0.75, hS),
                arrowprops=dict(arrowstyle="<->", color="#333", lw=1))
    ax.text(bS/2+0.80, hS/2, f"h={h*100:.0f}cm", fontsize=9, va="center", rotation=90)

    ax.set_xlim(-bS/2 - 0.9, bS/2 + 1.5)
    ax.set_ylim(-0.22, hS + 0.15)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Section — ferraillage réel", fontsize=11, fontweight="bold")

    plt.tight_layout()
    if nom_fichier:
        plt.savefig(nom_fichier, dpi=150, bbox_inches="tight")
    return fig


# ═══════════════════════════════════════════════════════════════
# 8.  DIAGRAMME DE DÉFORMATION — section + profil de déformation ELU
# ═══════════════════════════════════════════════════════════════

def diagramme_deformation(section, etat, nom_fichier=None):
    """
    Trace le schéma de section (avec ferraillage) et, à côté, le profil
    linéaire de déformation sur la hauteur à l'ELU (état renvoyé par
    etat_deformation_ELU), avec la position de l'axe neutre et les
    déformations aux nappes d'aciers.

    section : dict(b, h, c_inf, c_sup, As_inf, As_sup)
    etat    : dict renvoyé par etat_deformation_ELU() (non None)

    Retourne la figure matplotlib (fig).
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]
    As_inf, As_sup = section["As_inf"], section["As_sup"]

    fig, (ax_sec, ax_def) = plt.subplots(1, 2, figsize=(8, 6),
                                          gridspec_kw={"width_ratios": [1, 1.6]})

    # -- Section --
    sc = 1.0 / max(b, h)
    bS, hS = b * sc, h * sc
    ax_sec.add_patch(mpatches.Rectangle((-bS/2, -hS/2), bS, hS,
                                         fc="#e3e9f5", ec="#333", lw=1.5))
    r_b = 0.02
    y_sup = (h/2 - c_sup) * sc
    y_inf = (-h/2 + c_inf) * sc
    if As_sup > 0:
        ax_sec.add_patch(plt.Circle((0, y_sup), r_b, color="#b71c1c", zorder=4))
    if As_inf > 0:
        ax_sec.add_patch(plt.Circle((0, y_inf), r_b, color="#b71c1c", zorder=4))
    ax_sec.set_xlim(-0.6, 0.6) ;  ax_sec.set_ylim(-0.62, 0.62)
    ax_sec.set_aspect("equal") ;  ax_sec.axis("off")
    ax_sec.set_title("Section", fontsize=10, fontweight="bold")

    # -- Diagramme de déformation --
    y_top, y_bot = hS/2, -hS/2
    eps_top, eps_bot = etat["eps_top"], etat["eps_bot"]
    eps_max = max(abs(eps_top), abs(eps_bot), 1e-6)

    ax_def.axvline(0, color="#999", lw=0.8, ls=(0, (4, 3)))
    ax_def.plot([eps_top, eps_bot], [y_top, y_bot], color="#1565C0", lw=2.5)
    ax_def.plot([eps_top, eps_bot], [y_top, y_bot], "o", color="#1565C0", ms=6)

    # axe neutre : x est mesuré depuis la fibre comprimée (sup. si positif,
    # inf. sinon)
    if etat["positif"]:
        y_na = y_top - etat["x"] * sc
    else:
        y_na = y_bot + etat["x"] * sc
    ax_def.axhline(y_na, color="#2E7D32", lw=1, ls="--")
    ax_def.text(eps_max * 1.15, y_na, f"axe neutre\nx={etat['x']*1000:.0f}mm",
                fontsize=8, color="#2E7D32", va="center")

    ax_def.annotate(f"εsup={etat['eps_sup']*1e3:+.2f}‰", (eps_top, y_top),
                     xytext=(10, 8), textcoords="offset points", fontsize=8.5,
                     color="#1565C0", fontweight="bold")
    ax_def.annotate(f"εinf={etat['eps_inf']*1e3:+.2f}‰", (eps_bot, y_bot),
                     xytext=(10, -14), textcoords="offset points", fontsize=8.5,
                     color="#1565C0", fontweight="bold")

    ax_def.set_ylim(y_bot - 0.05, y_top + 0.05)
    ax_def.set_xlim(-eps_max * 1.4, eps_max * 1.4)
    ax_def.set_xlabel("Déformation ε", fontsize=9)
    ax_def.set_yticks([])
    ax_def.set_title("Diagramme de déformation (ELU, à la ruine)",
                      fontsize=10, fontweight="bold")
    ax_def.spines[["top", "right", "left"]].set_visible(False)

    plt.tight_layout()
    if nom_fichier:
        plt.savefig(nom_fichier, dpi=150, bbox_inches="tight")
    return fig


def lambda_eta_EC2(fck):
    """
    Coefficients λ (réduction de hauteur comprimée) et η (réduction de
    contrainte) du bloc rectangulaire équivalent, EC2 §3.1.7(3) :
      fck ≤ 50 MPa : λ=0,8  η=1,0
      fck > 50 MPa : λ=0,8-(fck-50)/400   η=1,0-(fck-50)/200
    """
    if fck <= 50.0:
        return 0.8, 1.0
    return 0.8 - (fck - 50.0) / 400.0, 1.0 - (fck - 50.0) / 200.0


def schema_bloc_rectangulaire(section, fck, fyk, etat, gamma_c=1.5, gamma_s=1.15,
                               nom_fichier=None):
    """
    Reproduit le schéma classique du cours (section / diagramme des
    déformations avec pivots A-B / bloc de contraintes rectangulaire
    équivalent Fc-Fs-z-Mu), rempli avec les valeurs réelles de l'état ELU
    calculé par etat_deformation_ELU().

    section : dict(b, h, c_inf, c_sup, As_inf, As_sup)
    etat    : dict renvoyé par etat_deformation_ELU() (non None)

    Retourne la figure matplotlib (fig).
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]
    As_inf, As_sup = section["As_inf"], section["As_sup"]

    mat = ec2.get_material_params(fck, fyk, gamma_c, gamma_s)
    fcd = mat["fcd"]
    lam, eta = lambda_eta_EC2(fck)

    positif = etat["positif"]
    d = (h - c_inf) if positif else (h - c_sup)
    x = etat["x"]
    eps_c = etat["eps_top"] if positif else etat["eps_bot"]
    eps_s = etat["eps_bot"] if positif else etat["eps_top"]

    lx = lam * x
    z = d - lx / 2.0
    Fc = eta * fcd * 1e3 * lx * b        # kN  (fcd MPa -> kN/m² = *1e3)
    As_t = As_inf if positif else As_sup
    fyd = mat["fyd"]
    # contrainte réelle dans l'acier tendu à l'état calculé (peut être
    # plastifiée ou non selon eps_s vs eps_yd)
    sigma_s = ec2.sigma_s(abs(eps_s), mat) if As_t > 0 else fyd
    Fs = As_t * sigma_s * 1e3            # kN
    Mu = Fc * z                          # kN.m (résultante béton, cohérence ELU)

    fig, axes = plt.subplots(1, 3, figsize=(12, 5),
                              gridspec_kw={"width_ratios": [0.8, 1, 1.3]})
    ax_sec, ax_def, ax_bloc = axes

    # ── Panneau 1 : Section ──────────────────────────────────────────────
    sc = 1.0 / h
    bS, hS = b * sc, h * sc
    ax_sec.add_patch(mpatches.Rectangle((-bS/2, 0), bS, hS, fc="white", ec="#333", lw=1.5))
    y_acier = c_inf * sc if positif else (h - c_sup) * sc
    ax_sec.add_patch(mpatches.Rectangle((-bS/2*0.55, y_acier - 0.012), bS*0.55, 0.024,
                                         fc="#334", ec="none"))
    ax_sec.annotate("", xy=(-bS/2-0.06, hS if positif else h*sc - y_acier),
                     xytext=(-bS/2-0.06, y_acier if positif else 0),
                     arrowprops=dict(arrowstyle="<->", color="#2E7D32", lw=1))
    ax_sec.text(-bS/2-0.10, (hS+y_acier)/2 if positif else y_acier/2, "d",
                fontsize=11, ha="right", va="center", style="italic")
    ax_sec.annotate("", xy=(-bS/2, -0.06), xytext=(bS/2, -0.06),
                     arrowprops=dict(arrowstyle="<->", color="#333", lw=1))
    ax_sec.text(0, -0.11, "b", ha="center", fontsize=11, style="italic")
    ax_sec.annotate("", xy=(bS/2+0.09, 0), xytext=(bS/2+0.09, hS),
                     arrowprops=dict(arrowstyle="<->", color="#333", lw=1))
    ax_sec.text(bS/2+0.13, hS/2, "h", fontsize=11, va="center", style="italic")
    ax_sec.set_xlim(-0.65, 0.75) ;  ax_sec.set_ylim(-0.18, hS+0.1)
    ax_sec.set_aspect("equal") ;  ax_sec.axis("off")

    # ── Panneau 2 : Diagramme des déformations (pivots A/B) ─────────────
    y_top, y_bot = hS, 0.0
    y_x = hS - x * sc if positif else x * sc  # position de l'axe neutre (échelle section)
    ax_def.axhline(y_x, color="#999", lw=0.8, ls=(0, (4, 3)))
    ax_def.annotate("", xy=(0, y_top), xytext=(0, y_x),
                     arrowprops=dict(arrowstyle="<->", color="#2E7D32", lw=1))
    ax_def.text(-0.10, (y_top+y_x)/2, "x", fontsize=11, ha="right", va="center",
                style="italic", color="#2E7D32")

    eps_max = max(abs(eps_c), abs(eps_s), 1e-6)
    e_axis_h = 0.9
    def e_to_x(eps):
        return eps / eps_max * e_axis_h

    ax_def.axvline(0, color="#666", lw=0.8)
    ax_def.plot([e_to_x(eps_c), e_to_x(eps_s)], [y_top, y_bot], color="#1565C0", lw=2.2)
    ax_def.plot(e_to_x(eps_c), y_top, "o", color="#1565C0", ms=6)
    ax_def.plot(e_to_x(eps_s), y_bot, "o", color="#1565C0", ms=6)
    ax_def.text(e_to_x(eps_c)+0.05, y_top, f"B\nεc={eps_c*1e3:.2f}‰",
                fontsize=9, color="#1565C0", va="center")
    ax_def.text(e_to_x(eps_s)+0.05, y_bot, f"A\nεs={eps_s*1e3:.2f}‰",
                fontsize=9, color="#1565C0", va="center")
    ax_def.annotate("", xy=(e_axis_h+0.15, hS/2), xytext=(-e_axis_h-0.15, hS/2),
                     arrowprops=dict(arrowstyle="->", color="#333", lw=1))
    ax_def.text(e_axis_h+0.20, hS/2, "ε", fontsize=11, style="italic", va="center")
    ax_def.annotate("", xy=(0, y_top+0.12), xytext=(0, y_bot-0.02),
                     arrowprops=dict(arrowstyle="->", color="#333", lw=1))
    ax_def.text(0.03, y_top+0.14, "y", fontsize=11, style="italic")
    ax_def.set_xlim(-e_axis_h-0.35, e_axis_h+0.4) ;  ax_def.set_ylim(-0.15, y_top+0.22)
    ax_def.axis("off")
    ax_def.set_title("Diagramme des déformations", fontsize=10, fontweight="bold")

    # ── Panneau 3 : Bloc rectangulaire équivalent (Fc, Fs, z, Mu) ────────
    lxS = lam * x * sc
    ax_bloc.add_patch(mpatches.Rectangle((0.05, hS - lxS), 0.35, lxS,
                                          fc="none", ec="#333", lw=1.3))
    ax_bloc.text(0.225, hS + 0.05, f"σc=η.fcd={eta*fcd:.1f} MPa",
                 fontsize=8.5, ha="center")
    ax_bloc.annotate("", xy=(0.42, hS), xytext=(0.42, hS - lxS),
                      arrowprops=dict(arrowstyle="<->", color="#2E7D32", lw=1))
    ax_bloc.text(0.46, hS - lxS/2, f"λ.x\n={lx*1000:.0f}mm", fontsize=8,
                 color="#2E7D32", va="center")

    y_Fc = hS - lxS/2
    ax_bloc.annotate("", xy=(0.05, y_Fc), xytext=(-0.15, y_Fc),
                      arrowprops=dict(arrowstyle="<-", color="#C62828", lw=2))
    ax_bloc.text(-0.18, y_Fc, f"Fc={Fc:.0f}kN", fontsize=9, color="#C62828",
                 ha="right", va="center", fontweight="bold")

    y_acier2 = c_inf * sc if positif else (h - c_sup) * sc
    ax_bloc.plot([0.05, 0.40], [y_acier2, y_acier2], color="#334", lw=3)
    ax_bloc.annotate("", xy=(0.55, y_acier2), xytext=(0.05, y_acier2),
                      arrowprops=dict(arrowstyle="->", color="#C62828", lw=2))
    ax_bloc.text(0.58, y_acier2, f"Fs={Fs:.0f}kN\n(σs={sigma_s:.0f}MPa)",
                 fontsize=8.5, color="#C62828", va="center", fontweight="bold")

    ax_bloc.annotate("", xy=(0.85, y_Fc), xytext=(0.85, y_acier2),
                      arrowprops=dict(arrowstyle="<->", color="#333", lw=1))
    ax_bloc.text(0.90, (y_Fc+y_acier2)/2, f"z={z*1000:.0f}mm", fontsize=9,
                 va="center", style="italic")

    ax_bloc.annotate("", xy=(0.20, y_acier2 - 0.08), xytext=(0.30, y_acier2 - 0.08),
                      arrowprops=dict(arrowstyle="->", color="#7B1FA2", lw=1.8,
                                       connectionstyle="arc3,rad=0.6"))
    ax_bloc.text(0.25, y_acier2 - 0.16, f"Mu={Mu:.0f} kN·m", fontsize=9.5,
                 color="#7B1FA2", ha="center", fontweight="bold")

    ax_bloc.set_xlim(-0.55, 1.15) ;  ax_bloc.set_ylim(-0.28, hS+0.15)
    ax_bloc.axis("off")
    ax_bloc.set_title("Bloc de contraintes équivalent", fontsize=10, fontweight="bold")

    fig.suptitle(f"Justification ELU flexion simple — pivot {'B' if positif else 'B (mirroir)'}"
                 f" — x={x*1000:.0f}mm, λ={lam:.3f}, η={eta:.3f}",
                 fontsize=10.5, fontweight="bold", y=1.02)
    plt.tight_layout()
    if nom_fichier:
        plt.savefig(nom_fichier, dpi=150, bbox_inches="tight")
    return fig


# ═══════════════════════════════════════════════════════════════
# 9.  RAPPORT PDF — entrants, sortants, figures, date
# ═══════════════════════════════════════════════════════════════

def generer_rapport_pdf(resultats, section, fck, fyk, M_ELS_kNm, gamma_c=1.5, gamma_s=1.15,
                         geom_inf=None, geom_sup=None,
                         nom_projet="", partie_ouvrage="",
                         nom_fichier="rapport_flexion_simple_EC2.pdf"):
    """
    Génère un rapport PDF multi-pages : page de garde (projet, partie
    d'ouvrage — optionnels, vides par défaut —, date, verdict), tableau
    des hypothèses (entrants), tableau des résultats (sortants), puis
    les figures (schéma de section détaillé, diagramme de déformation,
    bloc rectangulaire équivalent).

    geom_inf, geom_sup : dict renvoyés par geometrie_nappe() (celui qui
    n'est pas None doit correspondre à la nappe réellement tendue —
    passer None pour l'autre). Si aucun des deux n'est fourni, le schéma
    de section détaillé est omis du rapport.

    Retourne le chemin du fichier généré.
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import datetime

    if not nom_fichier.lower().endswith(".pdf"):
        nom_fichier += ".pdf"

    date_str = datetime.date.today().strftime("%d/%m/%Y")
    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]
    As_inf, As_sup = section["As_inf"], section["As_sup"]

    mr   = resultats["moment_reduit"]
    at_u = resultats["acier_theorique"]
    at_s = resultats["acier_theorique_ELS"]
    elu  = resultats["ELU"]
    els  = resultats["ELS_contraintes"]
    amin = resultats["aciers_minimaux"]
    fiss = resultats["fissuration"]

    def tag(ok):
        return "OK" if ok else "NON VERIFIE"

    with PdfPages(nom_fichier) as pdf:
        # ── Page 1 : garde + hypothèses + résultats ─────────────────────
        fig1 = plt.figure(figsize=(8.27, 11.69), facecolor="white")
        fig1.text(0.5, 0.965, "JUSTIFICATION À LA FLEXION SIMPLE",
                   fontsize=16, fontweight="bold", ha="center")
        fig1.text(0.5, 0.945, "Eurocode 2 (NF EN 1992-1-1 + Annexe Nationale française)",
                   fontsize=10, ha="center", color="#555")

        entete = (f"Projet            : {nom_projet if nom_projet else '—'}\n"
                  f"Partie d'ouvrage  : {partie_ouvrage if partie_ouvrage else '—'}\n"
                  f"Date              : {date_str}")
        fig1.text(0.08, 0.905, entete, fontsize=9.5, family="monospace", va="top")

        verdict_ok = resultats["verifie_global"]
        fig1.text(0.08, 0.83,
                  "✔  SECTION JUSTIFIÉE" if verdict_ok else "✘  SECTION NON JUSTIFIÉE",
                  fontsize=13, fontweight="bold",
                  color="#2E7D32" if verdict_ok else "#C62828")

        hyp = (
            f"BÉTON / ACIER\n"
            f"  Classe béton                 fck = {fck:.0f} MPa\n"
            f"  Acier                         fyk = {fyk:.0f} MPa\n"
            f"  Coefficients partiels         γc = {gamma_c:.2f}   γs = {gamma_s:.2f}\n"
            f"  Coeff. équivalence ELS        n = {els['detail']['n']:.2f}\n\n"
            f"GÉOMÉTRIE\n"
            f"  Section rectangulaire         b = {b*100:.0f} cm   h = {h*100:.0f} cm\n"
            f"  Enrobages                     c_inf = {c_inf*1000:.0f} mm   c_sup = {c_sup*1000:.0f} mm\n\n"
            f"FERRAILLAGE RÉEL\n"
            f"  As_inf mis en place           {As_inf*1e4:.2f} cm²\n"
            f"  As_sup mis en place           {As_sup*1e4:.2f} cm²\n\n"
            f"SOLLICITATIONS\n"
            f"  M_Ed ELU                      {elu['M_Ed']:.1f} kN·m\n"
            f"  M_Ed ELS (comb. caract.)       {M_ELS_kNm:.1f} kN·m\n"
            f"  Classe d'exposition            {fiss['classe_exposition']}"
        )
        fig1.text(0.08, 0.78, "ENTRANTS (hypothèses)", fontsize=11, fontweight="bold")
        fig1.text(0.08, 0.755, hyp, fontsize=8.7, family="monospace", va="top", linespacing=1.5)

        res_txt = (
            f"MOMENT RÉDUIT\n"
            f"  µ = {mr['mu']:.4f}    µlim = {mr['mu_lim']:.4f}\n"
            f"  Aciers comprimés nécessaires : {'OUI' if mr['besoin_aciers_comprimes'] else 'NON'}\n\n"
            f"ACIER THÉORIQUE (dimensionnement)\n"
            f"  As théorique ELU               {at_u['As_theorique']*1e4:.2f} cm²  (cas {at_u['cas']})\n"
            f"  As théorique ELS               {at_s['As_theorique_ELS']*1e4:.2f} cm²  "
            f"(critère {at_s['critere_gouvernant']})\n\n"
            f"ACIERS MINIMAUX RÉGLEMENTAIRES\n"
            f"  As,min ELU (§9.2.1.1)          {amin['As_min_ELU']*1e4:.2f} cm²   [{tag(amin['verifie_ELU'])}]\n"
            f"  As,min ELS (§7.3.2)            {amin['As_min_ELS']*1e4:.2f} cm²   [{tag(amin['verifie_ELS'])}]\n\n"
            f"ELU — CAPACITÉ\n"
            f"  M_Rd = {elu['M_Rd']:.1f} kN·m   taux = {elu['taux']*100:.1f} %   [{tag(elu['verifie'])}]\n\n"
            f"ELS — CONTRAINTES (§7.2)\n"
            f"  σc = {els['sigma_c']:.2f} MPa  (lim. {els['sigma_c_lim']:.2f})   [{tag(els['verifie_beton'])}]\n"
            f"  σs = {els['sigma_s']:.1f} MPa  (lim. {els['sigma_s_lim']:.1f})   [{tag(els['verifie_acier'])}]\n\n"
            f"OUVERTURE DE FISSURE (§7.3.4)\n"
            f"  wk = {fiss['wk']:.3f} mm   (wmax = {fiss['wmax']:.2f} mm)   [{tag(fiss['verifie'])}]"
        )
        fig1.text(0.08, 0.44, "SORTANTS (résultats)", fontsize=11, fontweight="bold")
        fig1.text(0.08, 0.415, res_txt, fontsize=8.7, family="monospace", va="top", linespacing=1.5)

        fig1.text(0.5, 0.02, "Outil interne — vérifier avant utilisation en note de calcul.",
                  fontsize=7.5, ha="center", color="#777", style="italic")
        pdf.savefig(fig1) ;  plt.close(fig1)

        # ── Page 2 : schéma de section détaillé ─────────────────────────
        if geom_inf is not None or geom_sup is not None:
            fig2 = schema_section_detaille(b, h, geom_inf=geom_inf, geom_sup=geom_sup)
            fig2.suptitle("Schéma de section — ferraillage réel", fontsize=12, fontweight="bold")
            pdf.savefig(fig2) ;  plt.close(fig2)

        # ── Page 3 : diagramme de déformation + bloc équivalent ─────────
        etat = resultats["deformation_ELU"]
        if etat is not None:
            fig3 = diagramme_deformation(section, etat)
            pdf.savefig(fig3) ;  plt.close(fig3)

            fig4 = schema_bloc_rectangulaire(section, fck, fyk, etat, gamma_c, gamma_s)
            pdf.savefig(fig4) ;  plt.close(fig4)

    print(f"  Rapport PDF sauvegardé : {nom_fichier}")
    return nom_fichier
