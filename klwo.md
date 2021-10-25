---
title: Kolmogorov lens in wave optics
author: Job Feldbrugge and Neil Turok
header-includes:
    <link rel = "icon" href = "figures/icon.png" type = "image/x-icon"> 
    <div class="topnav">
      <a href="https://jfeldbrugge.github.io/">Home</a>
      <a href="https://jfeldbrugge.github.io/Projects-and-Codes/">Projects and Code</a>
    </div>
---

<!--- <a href="">Feldbrugge, and Turok (2020)</a> ---> 

We consider lensing by a plasma whose electron density takes the form of a Gaussian random field with a Kolmogorov power spectrum: in the thin lens approximation, the lens is effectively two-dimensional and the power spectrum $P(k)\propto k^{-{5\over 3}}$. This is a good model for twinkling produced by a turbulent medium.

The following picture shows a realization of the lens:

<figure>
<img src='figures/figures_klwo/Kolmogorov_2.png' width=50% />
<figcaption> Fig. 1- A 2d Gaussian random lens with a Kolmogorov power spectrum. </figcaption>
</figure>

We calculate the intensity pattern produced by this lens, comprising an intricate caustic network with diffraction fringes:

<figure>
<table align='left' width=100% id="FIG">
<tr>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 800_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 800_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 1600_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 1600_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 1600_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 1600_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 3200_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 3200_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 3200_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 3200_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 6400_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 6400_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 6400_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 6400_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 12800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=1_ν = 12800_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 12800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Kolmogorov_2_f=2_ν = 12800_All.png' width=100% /></td>
</tr>
 </table>
 <figcaption> Fig. 2- The intensity pattern produced by lensing of the 2d Kolmogorov lens, for various lens strengths.</figcaption>
</figure>

We also consider a 2d Gaussian random field with a  power spectrum $P(k)\propto k^{-1}$:

<figure>
<img src='figures/figures_klwo/Inverse.png' width=50% />
<figcaption> Fig. 3- A realization of a random lens with power spectrum $P(k)\propto k^{-1}$. </figcaption>
</figure>

The intensity pattern for this case are shown below, again for various lens strengths.

<figure>
<table align='left' width=100% id="FIG">
<tr>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 800_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 800_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 1600_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 1600_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 1600_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 1600_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 3200_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 3200_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 3200_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 3200_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 6400_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 6400_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 6400_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 6400_All.png' width=100% /></td>
</tr>
<tr>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 12800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=1_ν = 12800_All.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 12800_Automatic.png' width=100% /></td>
<td><img src='figures/figures_klwo/Inverse_f=2_ν = 12800_All.png' width=100% /></td>
</tr>
</table>
 <figcaption> Fig. 4- The intensity pattern produced by a 2d plasma lens with power spectrum $P(k)\propto k^{-1}$, for various strengths of the lens.</figcaption>
</figure>
