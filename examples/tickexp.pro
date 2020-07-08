FUNCTION tickexp,axis, index, value, level; makes a string that is only the 10^EXP part of the value

if (value eq 0d0) then return, '0' 
valstr = strng(value,format = '(E10.2)')
expstr = strng(fix(strmid(valstr,2,/rev)))
totstr = '10!E'+expstr+'!N'

return,totstr
END
