load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
test = (/"alf", "alfx", "alff"/)
nt = dimsizes(test)
;ntrunc = 39
nproc = 5
missing = -9.99d33
; define your data dir
datadir = "./"
nlat = (ntrunc+1)*3/2
begin
  wks = gsn_open_wks("png","checksum_T"+ntrunc)
  gsn_define_colormap(wks,"amwg_blueyellowred")

  res = True
  res@cnFillOn = True
  res@cnFillMode = "RasterFill"
  res@cnLinesOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000/)*1.0d-16
;  res@cnLevels = (/5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 100000/)*1.0d-16
  res@gsnSpreadColors = True
  res@lbOrientation = "Vertical"
  res@tmYLLabelFontHeightF = 0.03
  res@gsnStringFontHeightF = 0.03
  res@gsnLeftString = "T"+ntrunc
  res@gsnMaximize = True

  setfileoption("bin","ReadByteOrder","BigEndian")
  wavenum = ispan(0,ntrunc,1)
  do j=0, nt-1
    res@tiMainString = test(j)
    do i=0, nproc-1
      fname = datadir+ test(j)+"_T"+sprinti("%0.5d",ntrunc)+"J" \
        +sprinti("%0.5d",nlat)+"_checksum_"+sprinti("%0.2d",i)+".dat"
      print(""+fname)
      if (i.eq.0) then
        ppw = fbindirread(fname,0,(/ntrunc+1,ntrunc+1/),"double")
        ppw@_FillValue = missing
        do m=1, ntrunc
          ppw(m,0:m-1) = missing
        end do
      else
        ppw = ppw + fbindirread(fname,0,(/ntrunc+1,ntrunc+1/),"double")
      end if
    end do
    ppw = fabs(1.0d0-(/transpose(ppw)/))
    ppw = mask(ppw, ppw.gt.1, False)
    res@gsnRightString = sprintf("%5.2e",avg(ppw))+" "+sprintf("%5.2e",max(ppw))
    plot = gsn_csm_contour(wks,ppw,res)
  end do
end
