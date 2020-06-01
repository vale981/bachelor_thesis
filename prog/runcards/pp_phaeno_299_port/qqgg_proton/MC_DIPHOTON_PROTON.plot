# BEGIN PLOT /MC_DIPHOTON_PROTON/.*
XTwosidedTicks=1
YTwosidedTicks=1
LegendAlign=r
RatioPlotErrorBandColor=yellow
RatioPlot=1
# END PLOT

# BEGIN HISTOGRAM /MC_DIPHOTON_PROTON/.*
ErrorBars=1
# END HISTOGRAM

# BEGIN PLOT /MC_DIPHOTON_PROTON/pT
Title=$\pT$ of leading photon
XLabel=$\pT$ [GeV]
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/pT_subl
Title=$\pT$ of sub-leading photon
XLabel=$\pT$ [GeV]
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/inv_m
Title=$m_{\gamma\gamma}$ of the two leading photons
XLabel=$m_{\gamma\gamma}$ [GeV]
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/eta
Title=$\eta$ of leading photon
XLabel=$\eta$
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/cos_theta
Title=$\cos\theta$ of leading photon
XLabel=$\cos\theta$
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/azimuthal_angle
Title=Azimuthal angle between the photons
XLabel=$\pi-\Delta\phi_{\gamma\gamma}$ [rad]
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\phi_{\gamma\gamma}$ [fb rad$^{-1}$]
LogY=0
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/o_angle
Title=scattering angle
XLabel=$\cos\theta^\ast$
LogY=0
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/o_angle_cs
Title=scattering angle in CS frame
XLabel=$\cos\theta^\ast_\text{CS}$
LogY=0
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/total_pT
Title=total $\pT$
XLabel=$\pT$
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/xs
Title=total cross section
YLabel=$\sigma$
LogY=0
YMin=42
YMax=48
# END PLOT
