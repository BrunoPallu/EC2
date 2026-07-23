"""
Application web Streamlit — Diagramme d'interaction N-M (Eurocode 2)
======================================================================
Interface web pour diagramme_interaction_EC2.py (moteur de calcul,
inchangé). Déployable gratuitement sur Streamlit Community Cloud.

Fichier principal à déclarer dans les paramètres Streamlit Cloud :
    app.py
"""

import matplotlib
matplotlib.use("Agg")  # backend sans interface graphique (obligatoire sur serveur)

import io
import tempfile
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import diagramme_interaction_EC2 as ec2

st.set_page_config(page_title="Diagramme d'interaction EC2", layout="wide")

# ── Classes de béton normalisées — EC2 Tableau 3.1 ──────────────────────────
CLASSES_BETON = [
    ("C12/15", 12.0), ("C16/20", 16.0), ("C20/25", 20.0),
    ("C25/30", 25.0), ("C30/37", 30.0), ("C35/45", 35.0),
    ("C40/50", 40.0), ("C45/55", 45.0), ("C50/60", 50.0),
    ("C55/67", 55.0), ("C60/75", 60.0), ("C70/85", 70.0),
    ("C80/95", 80.0), ("C90/105", 90.0),
]

st.title("Diagramme d'interaction N-M — Eurocode 2")
st.caption("Sections rectangulaire et circulaire — méthode des pivots A/B/C "
           "(EC2 §6.1, cours ENPC BAEP1 Séance 4)")

# ═══════════════════════════════════════════════════════════════════════
# BARRE LATÉRALE — saisie des hypothèses
# ═══════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.header("1. Section")
    forme = st.radio("Type de section", ["rect", "circ"],
                      format_func=lambda x: "Rectangulaire" if x == "rect" else "Circulaire")

    st.header("2. Matériaux")
    label_beton = st.selectbox("Classe béton", [c[0] for c in CLASSES_BETON], index=3)
    fck = dict(CLASSES_BETON)[label_beton]
    fyk = st.number_input("fyk [MPa]", value=500.0, step=25.0)
    col_g1, col_g2 = st.columns(2)
    gamma_c = col_g1.number_input("γc", value=1.5, step=0.05, format="%.2f")
    gamma_s = col_g2.number_input("γs", value=1.15, step=0.05, format="%.2f")

    st.header("3. Géométrie & ferraillage")
    if forme == "rect":
        col1, col2 = st.columns(2)
        b = col1.number_input("b [m]", value=1.00, step=0.05, min_value=0.05)
        h = col2.number_input("h [m]", value=0.40, step=0.05, min_value=0.05)
        c_sup = col1.number_input("Enrobage sup. [m]", value=0.055, step=0.005, format="%.3f")
        c_inf = col2.number_input("Enrobage inf. [m]", value=0.055, step=0.005, format="%.3f")

        st.markdown("**Nappe supérieure**")
        col1, col2 = st.columns(2)
        nb_sup = col1.number_input("Nb barres sup.", value=5, step=1, min_value=0)
        phi_sup = col2.number_input("Ø sup. [mm]", value=10.0, step=1.0)

        st.markdown("**Nappe inférieure**")
        col1, col2 = st.columns(2)
        nb_inf = col1.number_input("Nb barres inf.", value=10, step=1, min_value=0)
        phi_inf = col2.number_input("Ø inf. [mm]", value=10.0, step=1.0)

        As_sup = nb_sup * np.pi * (phi_sup / 1000 / 2) ** 2
        As_inf = nb_inf * np.pi * (phi_inf / 1000 / 2) ** 2
        section = dict(b=b, h=h, c_inf=c_inf, c_sup=c_sup,
                        As_inf=As_inf, As_sup=As_sup,
                        nb_sup=nb_sup, nb_inf=nb_inf)
    else:
        col1, col2 = st.columns(2)
        D = col1.number_input("D [m]", value=0.60, step=0.05, min_value=0.05)
        c_enr = col2.number_input("Enrobage [m]", value=0.070, step=0.005, format="%.3f")
        col1, col2 = st.columns(2)
        nb_barres = col1.number_input("Nb barres", value=12, step=1, min_value=3)
        phi_barre = col2.number_input("Ø barres [mm]", value=16.0, step=1.0)

        As_tot = nb_barres * np.pi * (phi_barre / 1000 / 2) ** 2
        section = dict(D=D, c_enr=c_enr, nb_barres=nb_barres, As_tot=As_tot)

    st.header("4. Sollicitations à vérifier (optionnel)")
    sol_txt = st.text_area("Une ligne par cas : N [kN], M [kN·m], label",
                            placeholder="2189, 842, ELU1\n3000, 800, ELU2", height=90)

# ═══════════════════════════════════════════════════════════════════════
# CALCUL (relancé automatiquement à chaque changement de paramètre)
# ═══════════════════════════════════════════════════════════════════════
sollicitations = []
for ligne in sol_txt.strip().splitlines():
    if not ligne.strip():
        continue
    parts = [p.strip() for p in ligne.split(",")]
    if len(parts) >= 2:
        try:
            Ned, Med = float(parts[0]), float(parts[1])
            lbl = parts[2] if len(parts) >= 3 else f"Cas{len(sollicitations)+1}"
            sollicitations.append((Ned, Med, lbl))
        except ValueError:
            pass

try:
    (N_arr, M_arr, mat, pts_cles, fibres, arma, Ac, H, zones
     ) = ec2.diagramme_interaction(
        section_type=forme, section_params=section,
        fck=fck, fyk=fyk, gamma_c=gamma_c, gamma_s=gamma_s,
        n_div=300, n_piv_A=150, n_piv_B=200, n_piv_C=100)
except Exception as e:
    st.error(f"Erreur de calcul : {type(e).__name__} — {e}")
    st.stop()

# ═══════════════════════════════════════════════════════════════════════
# AFFICHAGE
# ═══════════════════════════════════════════════════════════════════════
col_diag, col_sec = st.columns([2.4, 1])

with col_diag:
    fig_plotly = ec2.tracer_interactif(
        N_arr, M_arr, mat, pts_cles, forme, section,
        fibres, arma, Ac, H, zones=zones,
        sollicitations=sollicitations if sollicitations else None,
        nom_fichier=tempfile.NamedTemporaryFile(suffix=".html", delete=False).name)
    st.plotly_chart(fig_plotly, width="stretch")

with col_sec:
    fig_sec = plt.figure(figsize=(4.2, 8.0), facecolor="#f4f6f9")
    gs = GridSpec(2, 1, height_ratios=[1.15, 1.85], hspace=0.08, figure=fig_sec)
    ax_sec = fig_sec.add_subplot(gs[0])
    ax_info = fig_sec.add_subplot(gs[1])
    ax_info.set_xlim(0, 1) ; ax_info.set_ylim(0, 1) ; ax_info.axis("off")
    ec2.dessiner_section(ax_sec, ax_info, forme, section, mat, arma, Ac)
    st.pyplot(fig_sec, width="stretch")

st.subheader("Points caractéristiques")
df_pts = pd.DataFrame(
    [(k.strip(), f"{v:,.1f}".replace(",", " ")) for k, v in pts_cles.items()],
    columns=["Grandeur", "Valeur"])
st.dataframe(df_pts, hide_index=True, width="stretch")

if sollicitations:
    st.subheader("Vérification des sollicitations")
    from matplotlib.path import Path as MplPath
    path = MplPath(list(zip(N_arr, M_arr)))
    lignes = []
    for (Ned, Med, lbl) in sollicitations:
        ok = path.contains_point((Ned, Med))
        lignes.append((lbl, f"{Ned:.0f}", f"{Med:.0f}",
                        "✔ Vérifié" if ok else "✘ Dépassement"))
    st.dataframe(pd.DataFrame(lignes, columns=["Cas", "N_Ed [kN]", "M_Ed [kN·m]", "Statut"]),
                 hide_index=True, width="stretch")

# ═══════════════════════════════════════════════════════════════════════
# EXPORTS
# ═══════════════════════════════════════════════════════════════════════
st.subheader("Exports")
col1, col2, col3 = st.columns(3)

with col1:
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        ec2.tracer(N_arr, M_arr, mat, pts_cles, forme, section,
                   fibres, arma, Ac, H, zones=zones,
                   sollicitations=sollicitations if sollicitations else None,
                   nom_fichier=tmp.name)
        plt.close("all")
        pdf_bytes = open(tmp.name, "rb").read()
    st.download_button("📄 Rapport PDF", data=pdf_bytes,
                        file_name="rapport_diagramme_EC2.pdf", mime="application/pdf")

with col2:
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        ec2.exporter_csv(N_arr, M_arr, zones=zones, nom_fichier=tmp.name)
        csv_bytes = open(tmp.name, "rb").read()
    st.download_button("📊 Coordonnées CSV", data=csv_bytes,
                        file_name="diagramme_EC2.csv", mime="text/csv")

with col3:
    html_bytes = fig_plotly.to_html(include_plotlyjs="cdn").encode("utf-8")
    st.download_button("🌐 Diagramme HTML", data=html_bytes,
                        file_name="diagramme_EC2.html", mime="text/html")

st.caption("Outil interne — vérifier les résultats avant utilisation en note de calcul.")
