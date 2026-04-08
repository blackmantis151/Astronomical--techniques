pro make_grid, x, y, a, grid

;print,'---------------------------------'
;print,'Syntax: make_grid, x, y, a, grid'
;print,'x=5, y=10, a=0.5, and grid is the output'
;print,'---------------------------------'

grida=fltarr(((y-x)/a)+1)

int=0
for i=0, n_elements(grida)-1 do begin
grida(i)=int+x
int=int+a
endfor

grid=grida

end 
