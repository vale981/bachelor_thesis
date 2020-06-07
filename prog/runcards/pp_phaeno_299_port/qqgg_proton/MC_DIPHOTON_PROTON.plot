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
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\text{T}^{\gamma_1}$ [pb GeV$^{-1}$]
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/pT_subl
Title=$\pT$ of sub-leading photon
XLabel=$\pT$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\text{T}^{\gamma_2}$ [pb GeV$^{-1}$]
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/inv_m
Title=$m_{\gamma\gamma}$ of the two leading photons
XLabel=$m_{\gamma\gamma}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{\gamma\gamma}$ [pb GeV$^{-1}$]
LogX=1
XMin=1
XMax=2000
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/eta
Title=$\eta$ of leading photon
XLabel=$\eta$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta$ [pb]
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/cos_theta
Title=$\cos\theta$ of leading photon
XLabel=$\cos\theta$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\cos_{\theta}$ [pb]
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/azimuthal_angle
Title=Azimuthal angle between the photons
XLabel=$\pi-\Delta\phi_{\gamma\gamma}$ [rad]
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\phi_{\gamma\gamma}$ [pb rad$^{-1}$]
LogY=0
LogX=0
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/o_angle
Title=scattering angle
XLabel=$|\cos\theta^\ast|$
YLabel=$\mathrm{d}\sigma/\mathrm{d}|\cos\theta^\ast|$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/o_angle_cs
Title=scattering angle in CS frame
XLabel=$|\cos\theta^\ast_\text{CS}|$
YLabel=$\mathrm{d}\sigma/\mathrm{d}|\cos\theta^\ast_\text{CS}|$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/total_pT
Title=$\pT$ of the $\gamma\gamma$ system
XLabel=$\pT$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\pT_{\gamma\gamma}$ [pb GeV$^{-1}$]
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_DIPHOTON_PROTON/xs
Title=total cross section
YLabel=$\sigma$ [pb]
LogY=0
YMin=32
YMax=40
XLabel=
LogY=0
XMinorTickMarks=0
LegendAlign=r
RatioPlot=0
# END PLOT
