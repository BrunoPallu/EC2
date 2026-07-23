import streamlit as st

st.title("🏗️ Outils de calcul béton armé — Eurocode 2")
st.caption("Bureau d'études — outils internes de justification, en remplacement des feuilles Excel")

st.markdown("---")

col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("📐 Flexion composée")
    st.markdown(
        "Diagramme d'interaction N-M (méthode des pivots A/B/C), "
        "sections **rectangulaire** et **circulaire**.\n\n"
        "- Diagramme interactif\n"
        "- Schéma de section + ferraillage\n"
        "- Export PDF, CSV, HTML"
    )
    st.page_link("page_flexion_composee.py", label="Ouvrir l'outil", icon="➡️")

with col2:
    st.subheader("📏 Flexion simple")
    st.markdown(
        "Justification ELU + ELS, section **rectangulaire**.\n\n"
        "- Capacité ELU (réutilise le diagramme d'interaction)\n"
        "- Contraintes ELS, ouverture de fissure\n"
        "- Aciers minimaux (ELU + fissuration)"
    )
    st.page_link("page_flexion_simple.py", label="Ouvrir l'outil", icon="➡️")

with col3:
    st.subheader("✂️ Effort tranchant")
    st.markdown(
        "Justification de l'effort tranchant (bielles, armatures "
        "transversales).\n\n"
        "*À construire — pas encore disponible.*"
    )
    st.page_link("page_tranchant.py", label="En savoir plus", icon="➡️")

st.markdown("---")
st.caption(
    "Calculs réalisés conformément à la NF EN 1992-1-1 et son Annexe "
    "Nationale française. Outil interne — vérifier les résultats avant "
    "utilisation en note de calcul."
)
