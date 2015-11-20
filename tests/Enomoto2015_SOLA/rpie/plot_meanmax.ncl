; 39 1.686067e-15 1.876565e-14 2.047343e-16 7.765095e-16 3.989401e-15 4.019007e-14 4.292862e-16 6.661338e-16 2.112103e-33 1.723707e-32 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
test = (/"alfx mean","alfx max","alff mean","alff max","alfq mean","alfq max"/)
datadir = "./"
begin
  nv = 11
  nt = 11

  buf = asciiread(datadir+"rpie2_meanmax.txt",(/nt,nv/),"double")
  x = buf(:,0)

  wks = gsn_open_wks("eps","rpie2_meanmax")
  res = True
;  res@gsnMaximize = True
  res@gsnDraw = False
  res@gsnFrame= False
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.8
  res@vpXF      = 1.0
  res@trXMinF = x(0)
  res@trXMaxF = x(nt-1)
  res@trYMinF = -16
  res@trYMaxF =  -6
  res@trXLog = True
  res@xyLineThicknessF = 4.0
  res@xyLineColors = (/"red","red","blue","blue","black","black"/)
  res@xyDashPatterns = (/1,0,1,0,1,0,1,0/)
  res@tmXBLabelsOn = False
  res@tmXBOn = False
;  res@pmLegendDisplayMode = "Always"
;  res@pmLegendSide = "Top"
;  res@pmLegendParallelPosF = 0.14
;  res@pmLegendOrthogonalPosF = -0.4
;  res@pmLegendWidthF = 0.15
;  res@pmLegendHeightF = 0.15
;  res@lgLabelFontHeightF = 0.02
;  res@lgItemPlacement = "ExplicitPlacement"
;  res@lgItemPositions = (/0.0, 0.2, 0.4, 0.6, 0.8, 1.0/)
;  res@lgPerimOn = False
;  res@lgItemOrder = (/1,0/)

  err = new((/nv-1, nt/),"double")
  do i = 0, nv - 2
    err(i,:) = log10(mask(buf(:,i+1),buf(:,i+1).le.0,False))
  end do
  plot = new(2, "graphic")
;  res@lgItemOrder = (/2, 3, 0 ,1/)
  res@xyExplicitLegendLabels = test(0:3)
  res@tiYAxisString = "log10 relative precision"
  plot(0) = gsn_csm_xy(wks, x, err(0:3,:), res)

  delete(res@tmXBLabelsOn)
  delete(res@tmXBOn)
  res@tmXBMode = "Explicit"
  res@tmXBValues = x
  res@tmXBLabels = x
  res@tmXBLabelAngleF = 45
  res@tiXAxisString = "truncation wave number"
  res@tiYAxisString = "log10 identity error"
  err(8,:) = err(8,:) + 20
  err(9,:) = err(9,:) + 20
;  delete(res@lgItemOrder)
;  res@pmLegendOrthogonalPosF = -0.5
;  res@lgItemOrder = (/4, 5, 0, 1, 2, 3/)
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels = test(0:5)
  plot(1) = gsn_csm_xy(wks, x, err(4:9, :), res)

  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelBottom = 0.1
  pres@gsnPanelRight = 0.9
  pres@gsnPanelFigureStringsFontHeightF = 0.025
  pres@gsnPanelFigureStrings = (/"(a)", "(b)"/)
  pres@gsnPanelFigureStringsPerimOn = False
  pres@amJust = "TopLeft"
  gsn_panel(wks, plot, (/2,1/), pres)
  
end
