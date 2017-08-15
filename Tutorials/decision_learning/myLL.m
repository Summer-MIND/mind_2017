function LL = myLL(p1, p2)
    LL = (p1.*p2)./(p1.*p2+((1-p1).*(1-p2))) ;
end