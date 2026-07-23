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


def contraintes_ELS(section, fck, M_ELS_kNm, Es=200_000.0):
    """
    Contraintes élastiques en section fissurée sous moment de service
    M_ELS (combinaison caractéristique, en kN·m — signe : positif = fibre
    sup. comprimée, aciers inf. tendus, comme pour le diagramme d'inter-
    action).

    Retourne dict(x, I_fiss, sigma_c, sigma_st, sigma_sc, n, Ecm)
    """
    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]
    As_inf, As_sup = section["As_inf"], section["As_sup"]
    d = h - c_inf

    carac = caracteristiques_beton(fck)
    Ecm = carac["Ecm"]
    n = Es / Ecm

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
                       coeffs=None):
    """
    Calcule wk selon (7.8)-(7.14).

    nb_barres_tendues, phi_barre_mm : nombre et diamètre des barres de la
    nappe tendue sous M_ELS (utilisés pour φ et l'espacement).
    duree_chargement : "court_terme" (kt=0,6) ou "long_terme" (kt=0,4).
    """
    if coeffs is None:
        coeffs = COEFFS_RECOMMANDES

    b, h = section["b"], section["h"]
    c_inf, c_sup = section["c_inf"], section["c_sup"]

    carac = caracteristiques_beton(fck)
    fct_eff = carac["fctm"]

    els = contraintes_ELS(section, fck, M_ELS_kNm, Es=Es)
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
                              duree_chargement="long_terme"):
    """
    Justification complète d'une section rectangulaire en flexion simple :
    capacité ELU, contraintes ELS, aciers minimaux (ELU + fissuration),
    ouverture de fissure wk vs wmax(classe d'exposition).

    Retourne un dict structuré avec tous les résultats + verdicts booléens.
    """
    resultats = {}

    # --- ELU ---
    elu = capacite_ELU_flexion_simple(section, fck, fyk, gamma_c, gamma_s)
    M_Rd = elu["M_Rd_pos"] if M_ELU_kNm >= 0 else elu["M_Rd_neg"]
    resultats["ELU"] = dict(
        M_Ed=M_ELU_kNm, M_Rd=M_Rd,
        verifie=abs(M_ELU_kNm) <= abs(M_Rd),
        taux=abs(M_ELU_kNm) / abs(M_Rd) if M_Rd != 0 else np.inf,
        detail=elu)

    # --- ELS : contraintes ---
    els = contraintes_ELS(section, fck, M_ELS_kNm)
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

    # --- Aciers minimaux ---
    amin_elu = acier_min_ELU(section, fck, fyk)
    amin_fiss = acier_min_fissuration(section, fck, fyk)
    As_inf = section["As_inf"] ;  As_sup = section["As_sup"]
    As_tendu = As_inf if M_ELU_kNm >= 0 else As_sup
    resultats["aciers_minimaux"] = dict(
        As_tendu_reel=As_tendu,
        As_min_ELU=amin_elu["As_min"], verifie_ELU=As_tendu >= amin_elu["As_min"],
        As_min_fissuration=amin_fiss["As_min"], verifie_fissuration=As_tendu >= amin_fiss["As_min"],
        detail_ELU=amin_elu, detail_fissuration=amin_fiss)

    # --- Ouverture de fissure ---
    wmax = WMAX_TABLE.get(classe_exposition, 0.3)
    fiss = ouverture_fissure(section, fck, M_ELS_kNm, nb_barres_tendues,
                              phi_barre_mm, duree_chargement=duree_chargement)
    resultats["fissuration"] = dict(
        wk=fiss["wk"], wmax=wmax, verifie=fiss["wk"] <= wmax,
        classe_exposition=classe_exposition, detail=fiss)

    resultats["verifie_global"] = all([
        resultats["ELU"]["verifie"],
        resultats["ELS_contraintes"]["verifie_beton"],
        resultats["ELS_contraintes"]["verifie_acier"],
        resultats["aciers_minimaux"]["verifie_ELU"],
        resultats["aciers_minimaux"]["verifie_fissuration"],
        resultats["fissuration"]["verifie"],
    ])

    return resultats


def imprimer_rapport(resultats):
    """Affiche un rapport texte lisible du dict retourné par justifier_flexion_simple()."""
    def tag(ok):
        return "✔ VÉRIFIÉ    " if ok else "✘ NON VÉRIFIÉ"

    print("═" * 68)
    print("  JUSTIFICATION FLEXION SIMPLE — EC2")
    print("═" * 68)

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
    print(f"  As,min fissuration (§7.3.2)    = {a['As_min_fissuration']*1e4:6.2f} cm²   {tag(a['verifie_fissuration'])}")

    f = resultats["fissuration"]
    print(f"\n[Ouverture de fissure]  classe {f['classe_exposition']}")
    print(f"  wk = {f['wk']:.3f} mm  ≤  wmax = {f['wmax']:.2f} mm   {tag(f['verifie'])}")
    print(f"  (formule : {f['detail']['formule']})")

    print("\n" + "═" * 68)
    verdict = "✔  SECTION JUSTIFIÉE" if resultats["verifie_global"] else "✘  SECTION NON JUSTIFIÉE"
    print(f"  {verdict}")
    print("═" * 68)
