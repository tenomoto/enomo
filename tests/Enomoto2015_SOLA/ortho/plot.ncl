load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
test = (/"alfx", "alff", "alfq"/)
nt = dimsizes(test)
;ntrunc = 1279
ntrunc = 39
datadir = "./"
begin
  buf = asciiread(datadir+"ortho_T"+ntrunc+".txt",(/ntrunc,nt+1/),"double")
  x = buf(:,0)
  y = transpose(buf(:,1:))
  y = log10(mask(y, y.le.0, False))

  wks = gsn_open_wks("png","ortho_T"+ntrunc)

  res = True
  res@gsnMaximize = True
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.8
  res@vpXF      = 1.0
  res@trXMinF = x(0)
  res@trXMaxF = x(ntrunc-2)
  res@trYMinF = -16
  res@trYMaxF =  -13
  res@tiXAxisString = "order m"
  res@tiYAxisString = "log10 max orthogonality error"
;  res@tmYLMode = "Explicit"
;  res@tmYLValues = ispan(res@trYMinF, res@trYMaxF, 1)
;  res@tmXBLabels = res@tmYLValues
  res@xyMonoDashPattern = True
  res@xyLineColors = (/"red","blue","black"/)
  res@xyLineThicknessF = 4.0
;  res@pmLegendDisplayMode = "Always"
;  res@pmLegendSide = "Top"
;  res@pmLegendParallelPosF = 0.9
;  res@pmLegendOrthogonalPosF = -0.4
;  res@pmLegendWidthF = 0.15
;  res@pmLegendHeightF = 0.12
;  res@lgLabelFontHeightF = 0.02
;  res@lgItemOrder = (/2, 1,0/)
;  res@lgItemPlacement = "ExplicitPlacement"
;  res@lgItemPositions = (/0.5, 0.75, 1.0/)
;  res@lgPerimOn = False

;  res@xyExplicitLegendLabels = test
  plot = gsn_csm_xy(wks, x, y, res)

end
