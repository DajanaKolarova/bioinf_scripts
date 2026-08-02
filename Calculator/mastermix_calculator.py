from ctypes import c_buffer

import streamlit as st

st.title("Laboratory Calculator for HPLC Screening")
st.markdown(
    "This tool automatically calculates cascaded master mixes, nucleotide divisions by ratios, and accounts for pipetting error margins.")

# 1. INPUT PARAMETERS
st.header("1. Experiment Settings")
col1, col2, col3 = st.columns(3)

with col1:
    v_reaction = st.number_input("Reaction volume (uL)", value=50.0, step=1.0)
with col2:
    reserve = st.number_input("Pipetting error margin (%)", value=10.0, step=1.0)
with col3:
    c_buffer_stock = st.number_input("Buffer stock (e.g., 10x)", value=10, step=1)

reserve_factor = 1 + (reserve / 100.0)

# --- 2. ENTER REACTION COUNTS FOR NUCLEOTIDES ---
st.header("2. Enter reaction counts for individual substrates and conditions")
st.markdown("Enter the actual number of vials you want to pipette.")

# Define basic user inputs
substrate = ["ATP", "ADP", "AMP"]
conditions = ["No enzyme (Blank)", "Low cat. (100 nM)", "High cat. (500 nM)"]

# Create a dynamic table (data structure) for the user
uzivatelske_vstupy = {}
for s in substrate:
    st.subheader(f"substrate: {s}")
    uzivatelske_vstupy[s] = {}
    cols = st.columns(len(conditions))
    for i, podm in enumerate(conditions):
        with cols[i]:
            uzivatelske_vstupy[s][podm] = st.number_input(
                f"{podm}", min_value=0, value=3 if s == "ATP" and i > 0 else 2, key=f"{s}_{podm}"
            )

# --- 3. COMPUTATIONAL LOGIC ---
# Step A: Sum of actual reactions
reakce_atp = sum(uzivatelske_vstupy["ATP"].values())
reakce_adp = sum(uzivatelske_vstupy["ADP"].values())
reakce_amp = sum(uzivatelske_vstupy["AMP"].values())
celkem_reakci = reakce_atp + reakce_adp + reakce_amp

# Step B: Generation of the pipetting protocol
st.divider()
st.header("📋 Resulting Pipetting Protocol")

if celkem_reakci == 0:
    st.warning("No reactions have been entered yet.")
else:
    # 1. Common Master Mix Base (Buffer + Water + DTT)
    # Calculated on the total sum of all reactions + margin
    celkovo_s_rezervou = celkem_reakci * reserve_factor

    st.subheader("Step 1: Master Mix Common Base (Buffer, Water, DTT)")
    st.markdown(f"Total equivalent reactions (incl. {reserve}% margin): **{celkovo_s_rezervou:.1f}**")

    obj_pufr_1 = (v_reaction / c_pufr_stock) * celkovo_s_rezervou
    obj_voda_1 = (v_reaction * 0.6) * celkovo_s_rezervou  # Example ratio for illustration
    obj_dtt_1 = (v_reaction * 0.05) * celkovo_s_rezervou

    st.info(f"""
    Pipette into the common tube:
    * **Buffer ({c_pufr_stock}x):** {obj_pufr_1:.1f} uL
    * **Water (HPLC):** {obj_voda_1:.1f} uL
    * **DTT:** {obj_dtt_1:.1f} uL
    """)

    # 2. Branching into individual nucleotides
    st.subheader("Step 2: Branching into Sub-Master Mixes (Nucleotides)")
    st.markdown(
        "Split the common base according to the reaction ratio for the given nucleotide and add the specific nucleotide.")

    for nukl, reakce_dict in uzivatelske_vstupy.items():
        pocet_nukl = sum(reakce_dict.values())
        if pocet_nukl > 0:
            potreba_s_rezervou = pocet_nukl * reserve_factor
            # How much base to take from the previous step
            zaklad_aliquot = (obj_pufr_1 + obj_voda_1 + obj_dtt_1) * (pocet_nukl / celkem_reakci)
            # How much nucleotide stock to add (assuming 1 mM stock -> 200 uM target = 1/5 of reaction volume)
            nukleotid_stock = v_reaction * 0.2 * potreba_s_rezervou

            st.write(f"👉 **Subgroup {nukl}** (Target for {pocet_nukl} reactions + margin):")
            st.code(
                f"• Take {zaklad_aliquot:.1f} uL from the common base\n• Add {nukleotid_stock:.1f} uL of {nukl} stock (1 mM)")