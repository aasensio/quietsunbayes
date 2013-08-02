pro test
	a = dblarr(8,20000)
	openr,2,'test.extract'
	readu,2,a
	close,2
	
	a = a[*,10000:*]
		
	!p.multi = [0,4,2,0,1]
	for i = 0, 7 do begin
		cgplot,a[i,*]
		print, mean(a[i,*]), stddev(a[i,*])
	endfor
	!p.multi = 0
	
	stop
end