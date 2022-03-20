hold on %графики все в одном окне
isd = 10
rad = isd / 3

xhex=[0]
yhex=[0] 

[towers_X, towers_Y] = addTowersRec2(2, isd, rad)

[centerX,centerY] = getc(mas,mas2,rad)

%plot(ggg,gggg2,'r*')
[nt,nt2] = gg(ggg,gggg2,rad, towers_X, towers_Y)
plot(nt,nt2,'b.')

function [towers_coord_X,towers_coord_Y] = addTowersRec2(numTiers, isd, rad)
    towers_coord_X = [0]
    towers_coord_Y = [0]    
    x=[1 1/2 -1/2 -1 -1/2 1/2] * isd
    y=[0 sqrt(3)/2 sqrt(3)/2 0 -sqrt(3)/2 -sqrt(3)/2] * isd
    for newTire = 1:numTiers
        for c = 1:length(towers_coord_X)
            offsetX = towers_coord_X(c)
            offsetY = towers_coord_Y(c)
            new_towers_group_coords_X = x - offsetX
            new_towers_group_coords_Y = y - offsetY
            [new_towers_X, new_towers_Y] = check_unique_coords(towers_coord_X, towers_coord_Y, new_towers_group_coords_X, new_towers_group_coords_Y)
            towers_coord_X = horzcat(towers_coord_X, new_towers_X)
            towers_coord_Y = horzcat(towers_coord_Y, new_towers_Y)
        end
    end
    
    for c = 1:length(towers_coord_X)
        offsetX = towers_coord_X(c)
        offsetY = towers_coord_Y(c)

        pgon1 = nsidedpoly(6, 'Center', [offsetX + rad/2, offsetY + sqrt(3)/2 * rad], 'Sidelength', rad)
        pgon2 = nsidedpoly(6, 'Center', [offsetX + rad/2, offsetY - sqrt(3)/2 * rad], 'Sidelength', rad)
        pgon3 = nsidedpoly(6, 'Center', [offsetX - rad, offsetY], 'Sidelength', rad)
        plot([pgon1, pgon2, pgon3])
    end
    
    plot(towers_coord_X,towers_coord_Y,'bO')
end

function [mas_X, mas_Y ] = check_unique_coords(towers_coord_X, towers_coord_Y, new_towers_group_coords_X, new_towers_group_coords_Y)
    mas_X = []
    mas_Y = []
    for i = 1:length(new_towers_group_coords_X)
        flag = 0
        curr_check_dot_X = new_towers_group_coords_X(i)
        curr_check_dot_Y = new_towers_group_coords_Y(i)
        for j = 1:length(towers_coord_X)
            if ((curr_check_dot_X == towers_coord_X(j)) && (curr_check_dot_Y == towers_coord_Y(j)))
                flag = 1
            end
        end
        
        if (flag == 0)
            mas_X = horzcat(mas_X, curr_check_dot_X)
            mas_Y = horzcat(mas_Y, curr_check_dot_Y)
        end
    end
    
end

function [g,gg] = getc(xp,yp,rad)
    g=[]
    gg=[]
    for c = 1:length(xp)
        oX1=xp(c)-rad
        oX2=xp(c)-rad/cos(60)/2
        oX3=oX2
        oY1=yp(c)
        oY2=yp(c)+rad*sqrt(3)/2
        oY3=yp(c)-rad*sqrt(3)/2
        g = horzcat(g,oX1)
        g = horzcat(g,oX2)
        g = horzcat(g,oX3)
        gg = horzcat(gg,oY1)
        gg = horzcat(gg,oY2)
        gg = horzcat(gg,oY3)
    end
end

function [gx,gy] = gg(xc,yc,rad, towers_X, towers_Y)
    N = 100; %Number of users
    R = rad; %Radius of Hexagon
    gy=[]
    gx=[]
    for c = 1:length(xc)
        v_x = (R * cos((0:6)*pi/3))+xc(c);
        v_y = (R * sin((0:6)*pi/3))+yc(c);

        c_x = (R-rand(1, 3*N)*2*R)+xc(c);
        c_y = (R-rand(1, 3*N)*2*R)+yc(c);

        IN = inpolygon(c_x, c_y, v_x, v_y);

        c_x = c_x(IN);
        c_y = c_y(IN);

        idx = randperm(length(c_x));
        c_x = c_x(idx(1:N));
        c_y = c_y(idx(1:N));
        gx = horzcat(gx, c_x)
        gy = horzcat(gy, c_y)
    end
end
