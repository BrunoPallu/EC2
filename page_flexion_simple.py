import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import streamlit as st

import flexion_simple_EC2 as fs

CLASSES_BETON = [
    ("C12/15", 12.0), ("C16/20", 16.0), ("C20/25", 20.0),
    ("C25/30", 25.0), ("C30/37", 30.0), ("C35/45", 35.0),
    ("C40/50", 40.0), ("C45/55", 45.0), ("C50/60", 50.0),
    ("C55/67", 55.0), ("C60/75", 60.0), ("C70/85", 70.0),
    ("C80/95", 80.0), ("C90/105", 90.0),
]

CLASSES_EXPOSITION = ["X0", "XC1", "XC2", "XC3", "XC4",
                      "XD1", "XD2", "XD3", "XS1", "XS2", "XS3"]

DIAMETRES_HA = [5, 6, 8, 10, 12, 14, 16, 20, 25, 32, 40]

st.title("Justification à la flexion simple — Eurocode 2")
st.caption("Section rectangulaire — ELU, ELS, aciers minimaux, ouverture de fissure "
           "(NF EN 1992-1-1 §6.1, §7.2, §7.3, §9.2.1.1 + Annexe Nationale française)")

# ═══════════════════════════════════════════════════════════════════════
# BARRE LATÉRALE — saisie des hypothèses
# ═══════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.header("1. Matériaux")
    label_beton = st.selectbox("Classe béton", [c[0] for c in CLASSES_BETON], index=3)
    fck = dict(CLASSES_BETON)[label_beton]
    fyk = st.number_input("fyk [MPa]", value=500.0, step=25.0)
    col_g1, col_g2 = st.columns(2)
    gamma_c = col_g1.number_input("γc", value=1.5, step=0.05, format="%.2f")
    gamma_s = col_g2.number_input("γs", value=1.15, step=0.05, format="%.2f")

    n_auto = st.checkbox("Calculer n=Es/Ecm(fck) automatiquement", value=False,
                          help="Décoché (défaut) : utilise la valeur forfaitaire n ci-dessous "
                               "(usuellement 15). Coché : calcul précis à partir de Ecm(fck), "
                               "variable selon la classe de béton.")
    n_valeur = st.number_input("Coefficient d'équivalence n = Es/Ec", value=15.0, step=1.0,
                                disabled=n_auto,
                                help="Valeur forfaitaire usuelle : 15 (béton ordinaire).")
    n_impose = None if n_auto else n_valeur

    st.header("2. Géométrie")
    col1, col2 = st.columns(2)
    b = col1.number_input("b [m]", value=0.30, step=0.05, min_value=0.05)
    h = col2.number_input("h [m]", value=0.60, step=0.05, min_value=0.05)
    c_inf = col1.number_input("Enrobage inf. [m]", value=0.045, step=0.005, format="%.3f")
    c_sup = col2.number_input("Enrobage sup. [m]", value=0.045, step=0.005, format="%.3f")

    st.header("3. Ferraillage (nappe tendue)")
    st.caption("Sous le moment ELU saisi ci-dessous : nappe inf. si M≥0, nappe sup. si M<0")
    phi_barre = st.selectbox("Ø barres", DIAMETRES_HA, index=DIAMETRES_HA.index(16),
                              format_func=lambda d: f"HA{d}")
    nb_lits = st.radio("Nombre de lits", [1, 2], horizontal=True)
    col1, col2 = st.columns(2)
    nb_lit1 = col1.number_input("Nb barres — lit 1 (près du parement)", value=4, step=1, min_value=1)
    nb_lit2 = col2.number_input("Nb barres — lit 2", value=2, step=1, min_value=0,
                                 disabled=(nb_lits == 1))
    entraxe_lits = st.number_input(
        "Entraxe vertical entre lits [mm]", value=float(phi_barre) + 20.0, step=1.0,
        disabled=(nb_lits == 1),
        help="Distance axe à axe entre lits. Défaut = Ø + 20mm (écartement usuel).")
    nb_barres_par_lit = [int(nb_lit1), int(nb_lit2) if nb_lits == 2 else 0]

    st.header("3bis. Aciers comprimés (optionnel)")
    avec_ac = st.checkbox("Définir des aciers comprimés", value=False,
                           help="Décoché par défaut : pas d'aciers comprimés (section simplement armée).")
    if avec_ac:
        phi_barre_c = st.selectbox("Ø barres comprimées", DIAMETRES_HA,
                                    index=DIAMETRES_HA.index(12), format_func=lambda d: f"HA{d}",
                                    key="phi_c")
        nb_lits_c = st.radio("Nombre de lits (comprimé)", [1, 2], horizontal=True, key="lits_c")
        col1, col2 = st.columns(2)
        nb_lit1_c = col1.number_input("Nb barres — lit 1", value=2, step=1, min_value=1, key="l1c")
        nb_lit2_c = col2.number_input("Nb barres — lit 2", value=0, step=1, min_value=0,
                                       disabled=(nb_lits_c == 1), key="l2c")
        entraxe_lits_c = st.number_input(
            "Entraxe vertical entre lits [mm]", value=float(phi_barre_c) + 20.0, step=1.0,
            disabled=(nb_lits_c == 1), key="ec")
        nb_barres_par_lit_c = [int(nb_lit1_c), int(nb_lit2_c) if nb_lits_c == 2 else 0]
    else:
        phi_barre_c = 12.0
        nb_barres_par_lit_c = [0, 0]
        entraxe_lits_c = phi_barre_c + 20.0

    st.header("4. Sollicitations")
    M_ELU = st.number_input("M_Ed ELU [kN·m]", value=150.0, step=10.0,
                             help="Positif = fibre inf. tendue, négatif = fibre sup. tendue")
    M_ELS = st.number_input("M_Ed ELS [kN·m] (comb. caractéristique)", value=105.0, step=10.0)

    st.header("5. Durabilité")
    classe_exposition = st.selectbox("Classe d'exposition", CLASSES_EXPOSITION, index=1)
    wmax_defaut = fs.WMAX_TABLE.get(classe_exposition, 0.3)
    col1, col2 = st.columns(2)
    wmax_inf = col1.number_input("wk,max fibre inf. [mm]", value=wmax_defaut, step=0.05,
                                  format="%.2f", help="Préréglé selon la classe d'exposition, modifiable.")
    wmax_sup = col2.number_input("wk,max fibre sup. [mm]", value=wmax_defaut, step=0.05,
                                  format="%.2f", help="Préréglé selon la classe d'exposition, modifiable.")
    duree = st.radio("Durée de chargement (fissuration)",
                      ["long_terme", "court_terme"],
                      format_func=lambda x: "Longue durée (kt=0,4)" if x == "long_terme"
                      else "Courte durée (kt=0,6)")

    st.header("6. Rapport (optionnel)")
    nom_projet = st.text_input("Projet", value="")
    partie_ouvrage = st.text_input("Partie d'ouvrage", value="")

# ═══════════════════════════════════════════════════════════════════════
# CALCUL
# ═══════════════════════════════════════════════════════════════════════
positif = M_ELU >= 0
cote = "inf" if positif else "sup"
cote_comp = "sup" if positif else "inf"
enrobage_nominal = c_inf if positif else c_sup
enrobage_comp = c_sup if positif else c_inf

geom = fs.geometrie_nappe(h, enrobage_nominal, nb_barres_par_lit, phi_barre,
                           entraxe_vertical_mm=entraxe_lits, cote=cote)
As = geom["As_total"]

if avec_ac:
    geom_c = fs.geometrie_nappe(h, enrobage_comp, nb_barres_par_lit_c, phi_barre_c,
                                 entraxe_vertical_mm=entraxe_lits_c, cote=cote_comp)
    As_c = geom_c["As_total"]
    c_eff_comp = geom_c["c_eff"] if As_c > 0 else enrobage_comp
else:
    geom_c = None
    As_c = 0.0
    c_eff_comp = enrobage_comp

# Section "effective" transmise au moteur de calcul : c_eff reproduit
# exactement le bras de levier réel (centroïde des lits), quel que soit
# le nombre de lits — aucune autre fonction du moteur n'a besoin d'être
# modifiée pour ça.
if positif:
    section = dict(b=b, h=h, c_inf=geom["c_eff"], c_sup=c_eff_comp, As_inf=As, As_sup=As_c)
else:
    section = dict(b=b, h=h, c_inf=c_eff_comp, c_sup=geom["c_eff"], As_inf=As_c, As_sup=As)

wmax_actif = wmax_inf if positif else wmax_sup

st.sidebar.metric("As nappe tendue (réel)", f"{As*1e4:.2f} cm²")
st.sidebar.caption(f"Bras de levier réel d = {geom['d_eff']*1000:.1f} mm "
                    f"(enrobage équivalent = {geom['c_eff']*1000:.1f} mm)")
if avec_ac:
    st.sidebar.metric("As' nappe comprimée (réel)", f"{As_c*1e4:.2f} cm²")

try:
    res = fs.justifier_flexion_simple(
        section, fck=fck, fyk=fyk,
        M_ELU_kNm=M_ELU, M_ELS_kNm=M_ELS,
        nb_barres_tendues=nb_barres_par_lit[0], phi_barre_mm=phi_barre,
        classe_exposition=classe_exposition,
        gamma_c=gamma_c, gamma_s=gamma_s,
        duree_chargement=duree,
        n_impose=n_impose,
        wmax_override=wmax_actif)
except Exception as e:
    st.error(f"Erreur de calcul : {type(e).__name__} — {e}")
    st.stop()

# ═══════════════════════════════════════════════════════════════════════
# AFFICHAGE — verdict global
# ═══════════════════════════════════════════════════════════════════════
if res["verifie_global"]:
    st.success("✔ Section justifiée — tous les critères sont satisfaits", icon="✅")
else:
    st.error("✘ Section NON justifiée — au moins un critère n'est pas satisfait", icon="⚠️")

st.subheader("Sections d'acier requises")
a0 = res["aciers_minimaux"]
at = res["acier_theorique"]
at_els = res["acier_theorique_ELS"]
col_a1, col_a2, col_a3, col_a4, col_a5 = st.columns(5)
col_a1.metric("As théorique ELU", f"{at['As_theorique']*1e4:.2f} cm²")
col_a2.metric("As théorique ELS", f"{at_els['As_theorique_ELS']*1e4:.2f} cm²")
col_a3.metric("As,min ELU (§9.2.1.1)", f"{a0['As_min_ELU']*1e4:.2f} cm²")
col_a4.metric("As,min ELS (§7.3.2)", f"{a0['As_min_ELS']*1e4:.2f} cm²")
col_a5.metric("As réel mis en place", f"{a0['As_tendu_reel']*1e4:.2f} cm²",
              delta="✔ suffisant" if (a0["verifie_ELU"] and a0["verifie_ELS"]) else "✘ insuffisant",
              delta_color="normal" if (a0["verifie_ELU"] and a0["verifie_ELS"]) else "inverse")
if at["cas"] == "double":
    st.warning(f"Double armature ELU : As tendu = As1({at['As1']*1e4:.2f}) + "
               f"As2({at['As2']*1e4:.2f}) = {at['As_theorique']*1e4:.2f} cm² — "
               f"As comprimé théorique = {at['As_comprime_theorique']*1e4:.2f} cm²",
               icon="⚠️")
st.caption(f"As théorique = section nécessaire pour équilibrer exactement M_Ed (dimensionnement) "
           f"— à distinguer des minima réglementaires, qui s'appliquent quel que soit M_Ed. "
           f"Critère gouvernant l'As théorique ELS : {at_els['critere_gouvernant']}.")

st.subheader("Moment réduit µ — besoin d'aciers comprimés")
mr = res["moment_reduit"]
col_mu1, col_mu2, col_mu3 = st.columns(3)
col_mu1.metric("µ (moment réduit)", f"{mr['mu']:.4f}")
col_mu2.metric("µlim (pivot B pur)", f"{mr['mu_lim']:.4f}")
if mr["besoin_aciers_comprimes"]:
    col_mu3.error("Aciers comprimés nécessaires", icon="⚠️")
else:
    col_mu3.success("Section simplement armée OK", icon="✅")

col1, col2 = st.columns(2)

with col1:
    st.subheader("ELU — Capacité en flexion")
    e = res["ELU"]
    st.metric("Taux de sollicitation", f"{e['taux']*100:.1f} %",
               delta=f"M_Rd = {e['M_Rd']:.1f} kN·m", delta_color="off")
    if e["verifie"]:
        st.success(f"M_Ed = {e['M_Ed']:.1f} kN·m ≤ M_Rd = {e['M_Rd']:.1f} kN·m")
    else:
        st.error(f"M_Ed = {e['M_Ed']:.1f} kN·m > M_Rd = {e['M_Rd']:.1f} kN·m")

    st.subheader("ELS — Limitation des contraintes (§7.2)")
    s = res["ELS_contraintes"]
    st.write(f"σc = **{s['sigma_c']:.2f} MPa**  (limite {s['sigma_c_lim']:.2f} MPa)  "
             + ("✔" if s["verifie_beton"] else "✘"))
    st.write(f"σs = **{s['sigma_s']:.1f} MPa**  (limite {s['sigma_s_lim']:.1f} MPa)  "
             + ("✔" if s["verifie_acier"] else "✘"))

with col2:
    st.subheader("Aciers minimaux")
    a = res["aciers_minimaux"]
    st.write(f"As réel = **{a['As_tendu_reel']*1e4:.2f} cm²**")
    st.write(f"As,min ELU (§9.2.1.1) = {a['As_min_ELU']*1e4:.2f} cm²  "
             + ("✔" if a["verifie_ELU"] else "✘"))
    st.write(f"As,min ELS / fissuration (§7.3.2) = {a['As_min_ELS']*1e4:.2f} cm²  "
             + ("✔" if a["verifie_ELS"] else "✘"))

    st.subheader("Ouverture de fissure (§7.3.4)")
    f = res["fissuration"]
    face_txt = "fibre inf." if positif else "fibre sup."
    st.write(f"wk = **{f['wk']:.3f} mm**  (wmax = {f['wmax']:.2f} mm, {face_txt}, "
             f"classe {f['classe_exposition']})  " + ("✔" if f["verifie"] else "✘"))
    st.caption(f"Formule retenue : {f['detail']['formule']}")

    if not f["verifie"]:
        if st.button("🔁 Calculer le nb de barres minimal pour respecter wk,max", key="btn_iter"):
            with st.spinner("Calcul itératif en cours..."):
                iter_res = fs.nb_barres_pour_fissure(
                    section, fck, M_ELS, phi_barre, wmax_actif,
                    gamma_c=gamma_c, gamma_s=gamma_s, n_impose=n_impose,
                    duree_chargement=duree)
            if iter_res["converge"]:
                st.success(f"**{iter_res['nb_barres_requis']} barres HA{phi_barre}** "
                           f"(As = {iter_res['As_requis']*1e4:.2f} cm²) suffisent : "
                           f"wk = {iter_res['wk_obtenu']:.3f} mm ≤ {wmax_actif:.2f} mm.\n\n"
                           f"Ajustez le nombre de barres du lit 1 dans la barre latérale "
                           f"en conséquence.", icon="✅")
            else:
                st.error(f"Pas de convergence avec des HA{phi_barre} seuls (60 barres testées "
                         f"sans succès) — essayez un diamètre plus fin ou plusieurs lits.",
                         icon="⚠️")

st.markdown("---")
st.subheader("Schéma ELU — section, déformations et bloc de contraintes (alignés)")
geom_inf_plot = geom if positif else geom_c
geom_sup_plot = geom_c if positif else geom

if res["deformation_ELU"] is not None:
    fig_complet = fs.schema_ELU_complet(section, fck, fyk, res["deformation_ELU"],
                                         gamma_c=gamma_c, gamma_s=gamma_s,
                                         geom_inf=geom_inf_plot, geom_sup=geom_sup_plot)
    st.pyplot(fig_complet, width="stretch")
    d = res["deformation_ELU"]
    col_i1, col_i2, col_i3 = st.columns(3)
    col_i1.metric("Axe neutre x", f"{d['x']*1000:.1f} mm")
    col_i2.metric("εsup", f"{d['eps_sup']*1e3:+.2f} ‰")
    col_i3.metric("εinf", f"{d['eps_inf']*1e3:+.2f} ‰")
    st.caption("Les 3 panneaux partagent la même échelle verticale : fibre supérieure, fibre "
               "inférieure et barycentre réel des aciers tendus sont à la même hauteur sur "
               "les 3 panneaux (repères horizontaux pointillés). État affiché à la ruine "
               "(pivot B, M=M_Rd) — pas l'état sous M_Ed de service. Le bras de levier tient "
               "compte du diamètre des barres, du nombre de lits et de l'enrobage réel.")
else:
    st.info("Pivot A probable (section très peu armée) — schéma non calculé pour ce cas ; "
            "affichage de la section seule.", icon="ℹ️")
    fig_sec = fs.schema_section_detaille(b, h, geom_inf=geom_inf_plot, geom_sup=geom_sup_plot)
    st.pyplot(fig_sec, width="stretch")

st.markdown("---")
st.subheader("Détail des résultats")

df = pd.DataFrame([
    ("µ (moment réduit)", f"{mr['mu']:.4f}"),
    ("µlim", f"{mr['mu_lim']:.4f}"),
    ("M_Rd (ELU)", f"{res['ELU']['M_Rd']:.1f} kN·m"),
    ("σc (ELS)", f"{res['ELS_contraintes']['sigma_c']:.2f} MPa"),
    ("σs (ELS)", f"{res['ELS_contraintes']['sigma_s']:.1f} MPa"),
    ("Axe neutre élastique (ELS) x", f"{res['ELS_contraintes']['detail']['x']*1000:.1f} mm"),
    ("As,min ELU", f"{res['aciers_minimaux']['As_min_ELU']*1e4:.2f} cm²"),
    ("As,min ELS / fissuration", f"{res['aciers_minimaux']['As_min_ELS']*1e4:.2f} cm²"),
    ("εsm (déformation moyenne acier)", f"{res['fissuration']['detail']['eps_sm']*1e3:.3f} ‰"),
    ("εcm (déformation moyenne béton tendu)", f"{res['fissuration']['detail']['eps_cm']*1e3:.3f} ‰"),
    ("εsm − εcm", f"{res['fissuration']['detail']['esm_ecm']*1e3:.3f} ‰"),
    ("wk", f"{res['fissuration']['wk']:.3f} mm"),
    ("sr,max", f"{res['fissuration']['detail']['sr_max']:.1f} mm"),
], columns=["Grandeur", "Valeur"])
st.dataframe(df, hide_index=True, width="stretch")

st.markdown("---")
st.subheader("Export")
import tempfile
geom_inf_pdf = geom if positif else geom_c
geom_sup_pdf = geom_c if positif else geom
with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
    fs.generer_rapport_pdf(
        res, section, fck, fyk, M_ELS, gamma_c=gamma_c, gamma_s=gamma_s,
        geom_inf=geom_inf_pdf, geom_sup=geom_sup_pdf,
        nom_projet=nom_projet, partie_ouvrage=partie_ouvrage,
        nom_fichier=tmp.name)
    pdf_bytes = open(tmp.name, "rb").read()
st.download_button("📄 Télécharger le rapport PDF", data=pdf_bytes,
                    file_name="rapport_flexion_simple_EC2.pdf", mime="application/pdf")

st.caption("Outil interne — vérifier les résultats avant utilisation en note de calcul. "
           "Valeurs conformes à la NF EN 1992-1-1/NA (mars 2007).")
