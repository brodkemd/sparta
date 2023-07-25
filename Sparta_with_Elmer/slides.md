---
theme: default
style: |
  section {
    background-color: #FFFFFF;
    font-color: black;
    font-size: 18pt;
  }
  h1 {
    color: #E00122;
    text-align: center;
    font-size: 30pt;
  }

---

<div style="color: #E00122; text-align: center; font-size: 60pt;"> Sparta + Elmer</div>
<div align="center"> <h3>Marek Brodke and Khare Prashant</h3></div>

---

# What each do Separately

__Sparta__:
- Simulates particles, making it excellent for rarefied flow.
- Simulate particle interactions with a solid surface.
- Can find forces and energy fluxes across a surface and its elements.

__Elmer__:
- Simulates solid body mechanics (i.e. deformation and heating).
- Can also simulate continuum flows, electro-dynamics, radiation, and more.

---

# What they do Together
By coupling Sparta with Elmer it is now possible to:
- Use particle interactions with a surface to find heating and deformation.
- Dynamically update surface mesh in both Sparta and Elmer to accommodate the results of heating and deformation.

---
# Diagram of the Operation

<div align="center">
<img src="Diagram.png" style="width:80%;height:auto;"></img>
</div>
