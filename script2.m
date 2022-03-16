hold on %графики все в одном окне
isd = 10
rad = isd / 3

%plot(pgon1)

xhex=[0]
yhex=[0] 
%pgon1 = nsidedpoly(6, 'Center', [rad/2, sqrt(3)/2 * rad], 'Sidelength', rad)

%pgon2 = nsidedpoly(6, 'Center', [rad/2, -sqrt(3)/2 * rad], 'Sidelength', rad)
%pgon3 = nsidedpoly(6, 'Center', [-rad, 0], 'Sidelength', rad)
%pg = plot(pgon1)
%pg.FaceColor = 'none';



mas = addTowersRec(xhex, yhex, 2, 10, rad)
disp('xplt');
disp(mas);
%pgon1 = nsidedpoly(6, 'Center', [rad/2, sqrt(3)/2 * rad], 'Sidelength', rad)

function  mas = addTowersRec(xplt, yplt, numTiers, isd, rad)
    x=[1 1/2 -1/2 -1 -1/2 1/2] * isd
    y=[0 sqrt(3)/2 sqrt(3)/2 0 -sqrt(3)/2 -sqrt(3)/2] * isd
    mas = []
    if (numTiers == 0)
        for c = 1:length(xplt)
            offsetX = xplt(c)
            offsetY = yplt(c)
            pgon1 = nsidedpoly(6, 'Center', [offsetX + rad/2, offsetY + sqrt(3)/2 * rad], 'Sidelength', rad)
            pgon2 = nsidedpoly(6, 'Center', [offsetX + rad/2, offsetY - sqrt(3)/2 * rad], 'Sidelength', rad)
            pgon3 = nsidedpoly(6, 'Center', [offsetX - rad, offsetY], 'Sidelength', rad)
            
            plot([pgon1, pgon2, pgon3])
        end
        plot(xplt,yplt,'bo')
    elseif (numTiers > 0)
        for c = 1:length(xplt)
            offsetX = xplt(c)
            offsetY = yplt(c)
            %mas(end+1) = offsetX
            xplt = horzcat(xplt, x-offsetX)
            yplt =  horzcat(yplt, y-offsetY)
            
            pgon1 = nsidedpoly(6, 'Center', [offsetX + rad/2, offsetY + sqrt(3)/2 * rad], 'Sidelength', rad)
            pgon2 = nsidedpoly(6, 'Center', [offsetX + rad/2, offsetY - sqrt(3)/2 * rad], 'Sidelength', rad)
            pgon3 = nsidedpoly(6, 'Center', [offsetX - rad, offsetY], 'Sidelength', rad)
            
            plot([pgon1, pgon2, pgon3])
            
            
            addTowersRec(xplt, yplt, numTiers-1, isd, rad)
        end
    end
end
