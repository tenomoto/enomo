load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
test = (/"alfx","alff", "alfx", "alff", "alfq"/)
datadir = "./"
begin
;  ntrunc = 10239 
  nlat = (ntrunc+1)*3/2
  nt = dimsizes(test)

  wks = gsn_open_wks("png","rpie2_T"+ntrunc)
  res = True
  res@gsnMaximize = True
  res@gsnDraw = False
  res@gsnFrame= False
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.8
  res@vpXF      = 0.1
  res@trXMinF =  0.
  res@trXMaxF = 90.
  res@trYMinF = -17
  res@trYMaxF =  -7
  res@xyMonoDashPattern = True
  res@xyLineColors = (/"red","blue","black"/)
  res@xyLineThicknessF = 4.0
  res@tmXBLabelsOn = False
  res@tmXBOn = False
;  res@pmLegendDisplayMode = "Always"
;  res@pmLegendSide = "Top"
;  res@pmLegendParallelPosF = 0.15
;  res@pmLegendOrthogonalPosF = -0.55
;  res@pmLegendWidthF = 0.15
;  res@pmLegendHeightF = 0.18
;  res@lgLabelFontHeightF = 0.02
;  res@lgItemPlacement = "ExplicitPlacement"
;  res@lgItemPositions = (/0.6, 0.8, 1.0/)
;  res@lgPerimOn = False
;  res@lgItemOrder = (/1,0/)

  i = 0
  err = new((/nt, nlat/2/),"double")
  buf = asciiread(datadir+"rpie2_T"+ntrunc+".txt",(/nlat/2,7/),"double")
  lat = buf(:,1)
  do i = 0, nt - 1
    err(i,:) = log10(mask(buf(:,i+2),buf(:,i+2).le.0,False))
  end do
;  print("rp: alfx="+max(err(0,:))+" alff="+max(err(1,:)))
;  print("ie: alfx="+max(err(2,:))+" alff="+max(err(3,:))+" alfq="+ max(err(4,:)))
  plot = new(2, "graphic")
;  res@lgItemOrder = (/1,0/)
  res@tiYAxisString = "log10 relative precision"
  res@xyExplicitLegendLabels = test(0:1)
;  res@trYMinF = -14
;  res@trYMaxF =  -5
  plot(0) = gsn_csm_xy(wks, lat, err(0:1,:), res)

  delete(res@tmXBLabelsOn)
  delete(res@tmXBOn)
  res@tmXBMode = "Manual"
  res@tmXBTickSpacingF = 15.
  res@tiXAxisString = "latitude"
  res@tiYAxisString = "log10 identity error"
;  res@trYMinF = -17
;  res@trYMaxF =  -7
  err(4,:) = err(4,:) + 20
;  delete(res@lgItemOrder)
;  res@lgItemOrder = (/2, 1,0/)
;  delete(res@xyExplicitLegendLabels)
;  res@xyExplicitLegendLabels = test(2:4)
  plot(1) = gsn_csm_xy(wks, lat, err(2:4,:), res)

  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelBottom = 0.15
  pres@gsnPanelFigureStrings = (/"(a)", "(b)"/)
  pres@gsnPanelFigureStringsPerimOn = False
  pres@gsnPanelFigureStringsFontHeightF = 0.025
  pres@amJust = "TopLeft"
  gsn_panel(wks, plot, (/2,1/), pres)
  
end
