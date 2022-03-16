isd = 10
n = 1
xhex=[0]
yhex=[0] 

disp('logs')
addTowersRec(xhex, yhex, 2, 10)



function  addTowersRec(xplt, yplt, numTiers, isd)
    x=[1 1/2 -1/2 -1 -1/2 1/2] * isd
    y=[0 sqrt(3)/2 sqrt(3)/2 0 -sqrt(3)/2 -sqrt(3)/2] * isd
    if (numTiers == 0)
        plot(xplt,yplt,'ro')
    elseif (numTiers > 0)
        for c = 1:length(xplt)
            offsetX = xplt(c)
            offsetY = yplt(c)
            xplt = horzcat(xplt, x-offsetX)
            yplt =  horzcat(yplt, y-offsetY)
            addTowersRec(xplt, yplt, numTiers-1, isd)
        end
    end
end
