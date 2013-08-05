pro test, n
	a = dblarr(408,n)
	openr,2,'test.extract'
	readu,2,a
	close,2
	
	a = a[*,n/2:*]
	
	stop
		
	!p.multi = [0,4,2,0,1]
	for i = 0, 7 do begin
		cgplot,a[i,*]
		print, mean(a[i,*]), stddev(a[i,*])
	endfor
	!p.multi = 0
	
	stop
end
