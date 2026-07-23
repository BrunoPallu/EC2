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

    st.header("2. Géométrie")
    col1, col2 = st.columns(2)
    b = col1.number_input("b [m]", value=0.30, step=0.05, min_value=0.05)
    h = col2.number_input("h [m]", value=0.60, step=0.05, min_value=0.05)
    c_inf = col1.number_input("Enrobage inf. [m]", value=0.045, step=0.005, format="%.3f")
    c_sup = col2.number_input("Enrobage sup. [m]", value=0.045, step=0.005, format="%.3f")

    st.header("3. Ferraillage (nappe tendue)")
    st.caption("Sous le moment ELU saisi ci-dessous : nappe inf. si M≥0, nappe sup. si M<0")
    col1, col2 = st.columns(2)
    nb_barres = col1.number_input("Nb barres", value=4, step=1, min_value=1)
    phi_barre = col2.number_input("Ø [mm]", value=16.0, step=1.0)

    st.header("4. Sollicitations")
    M_ELU = st.number_input("M_Ed ELU [kN·m]", value=150.0, step=10.0,
                             help="Positif = fibre inf. tendue, négatif = fibre sup. tendue")
    M_ELS = st.number_input("M_Ed ELS [kN·m] (comb. caractéristique)", value=105.0, step=10.0)

    st.header("5. Durabilité")
    classe_exposition = st.selectbox("Classe d'exposition", CLASSES_EXPOSITION, index=1)
    duree = st.radio("Durée de chargement (fissuration)",
                      ["long_terme", "court_terme"],
                      format_func=lambda x: "Longue durée (kt=0,4)" if x == "long_terme"
                      else "Courte durée (kt=0,6)")

# ═══════════════════════════════════════════════════════════════════════
# CALCUL
# ═══════════════════════════════════════════════════════════════════════
As = nb_barres * np.pi * (phi_barre / 1000 / 2) ** 2
if M_ELU >= 0:
    section = dict(b=b, h=h, c_inf=c_inf, c_sup=c_sup, As_inf=As, As_sup=0.0)
else:
    section = dict(b=b, h=h, c_inf=c_inf, c_sup=c_sup, As_inf=0.0, As_sup=As)

st.sidebar.metric("As nappe tendue", f"{As*1e4:.2f} cm²")

try:
    res = fs.justifier_flexion_simple(
        section, fck=fck, fyk=fyk,
        M_ELU_kNm=M_ELU, M_ELS_kNm=M_ELS,
        nb_barres_tendues=nb_barres, phi_barre_mm=phi_barre,
        classe_exposition=classe_exposition,
        gamma_c=gamma_c, gamma_s=gamma_s,
        duree_chargement=duree)
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
    st.write(f"wk = **{f['wk']:.3f} mm**  (wmax = {f['wmax']:.2f} mm, classe {f['classe_exposition']})  "
             + ("✔" if f["verifie"] else "✘"))
    st.caption(f"Formule retenue : {f['detail']['formule']}")

st.markdown("---")
st.subheader("Diagramme de déformation à l'ELU")
if res["deformation_ELU"] is not None:
    col_diag, col_txt = st.columns([1.6, 1])
    with col_diag:
        fig = fs.diagramme_deformation(section, res["deformation_ELU"])
        st.pyplot(fig, width="stretch")
    with col_txt:
        d = res["deformation_ELU"]
        st.write(f"Axe neutre x = **{d['x']*1000:.1f} mm**")
        st.write(f"εsup = **{d['eps_sup']*1e3:+.2f} ‰**")
        st.write(f"εinf = **{d['eps_inf']*1e3:+.2f} ‰**")
        st.caption("État à la ruine (pivot B), pas l'état sous M_Ed de service.")

    st.subheader("Schéma classique — bloc de contraintes équivalent")
    fig_bloc = fs.schema_bloc_rectangulaire(section, fck, fyk, res["deformation_ELU"],
                                             gamma_c=gamma_c, gamma_s=gamma_s)
    st.pyplot(fig_bloc, width="stretch")
    st.caption("Diagramme des déformations (pivots A/B) et bloc rectangulaire équivalent "
               "(λ, η — EC2 §3.1.7(3)), à titre de vérification croisée avec le calcul exact "
               "parabole-rectangle ci-dessus (Mu du bloc ≈ M_Rd exact à moins de 1% près).")
else:
    st.info("Pivot A probable (section très peu armée) — diagramme non calculé pour ce cas.",
            icon="ℹ️")

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
    ("wk", f"{res['fissuration']['wk']:.3f} mm"),
    ("sr,max", f"{res['fissuration']['detail']['sr_max']:.1f} mm"),
], columns=["Grandeur", "Valeur"])
st.dataframe(df, hide_index=True, width="stretch")

st.caption("Outil interne — vérifier les résultats avant utilisation en note de calcul. "
           "Valeurs conformes à la NF EN 1992-1-1/NA (mars 2007).")
