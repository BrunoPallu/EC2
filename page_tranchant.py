import streamlit as st

st.title("Effort tranchant — à venir")
st.info(
    "Cet outil n'est pas encore construit — on le fera ensemble dans une "
    "prochaine session.",
    icon="🚧",
)

st.markdown("""
D'après le cahier des charges, cet outil couvrira :

- Résistance à l'effort tranchant sans armatures (V_Rd,c)
- Dimensionnement des armatures d'effort tranchant (méthode des bielles,
  cotan(θ) variable, EC2 §6.2.3)
- Sections d'acier transversal minimales (ρw,min)
- Dispositions constructives (espacement des cadres/étriers)

Retour à l'[accueil](/page_accueil) une fois prêt à démarrer cet outil.
""")
