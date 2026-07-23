"""
Outils de calcul béton armé EC2 — Point d'entrée Streamlit
=============================================================
Application multi-pages : une page d'accueil qui distribue vers les
différents outils de justification (flexion composée, flexion simple,
effort tranchant à venir).

Fichier principal à déclarer dans les paramètres Streamlit Cloud :
    app.py
"""

import streamlit as st

st.set_page_config(page_title="Outils EC2", layout="wide", page_icon="🏗️")

accueil = st.Page("page_accueil.py", title="Accueil", icon="🏠", default=True)
flexion_composee = st.Page("page_flexion_composee.py", title="Flexion composée", icon="📐")
flexion_simple = st.Page("page_flexion_simple.py", title="Flexion simple", icon="📏")
tranchant = st.Page("page_tranchant.py", title="Effort tranchant", icon="✂️")

pg = st.navigation({
    "": [accueil],
    "Outils de justification": [flexion_composee, flexion_simple, tranchant],
})
pg.run()
