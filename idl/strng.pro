function strng,s,format=format
if not keyword_set(format) then format=''
return, strtrim(string(s,format=format),2)
end
