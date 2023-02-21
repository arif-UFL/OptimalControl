function dydt = myode1(t,y)
dydt = zeros(2,1);
dydt(1) = -y(1)-y(2);
dydt(2) = -y(1)+y(2);

end

    
    