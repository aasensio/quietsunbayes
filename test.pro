pro test, n
	a = dblarr(408,n)
	openr,2,'test.extract'
	readu,2,a
	close,2
	
; 	a = a[*,n/2:*]
	
	!p.multi = [0,8,8,0,1]
	for i = 0, 7 do begin
		for j = 0, 7 do begin
			if (i eq j) then begin
				cgplot,a[i,*]
				print, mean(a[i,*]), stddev(a[i,*])
			endif else begin
				cgplot, a[i,*], a[j,*], psym=3
				cgoplot, a[i,n/2:*], a[j,n/2:*], psym=3, col='red'
			endelse
		endfor
	endfor
	!p.multi = 0
	
	stop
end
