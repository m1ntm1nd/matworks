hold on %графики все в одном окне
isd = 10
rad = isd / 3

xhex=[0]
yhex=[0] 

[towers_X, towers_Y] = addTowersRec2(1, isd, rad)



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
    
    plot(towers_coord_X,towers_coord_Y, 'Or')

    users_X = []
    users_Y = []

    
    for c = 1:length(towers_coord_X)
        offsetX = towers_coord_X(c)
        offsetY = towers_coord_Y(c)
        center_1_X = offsetX + rad/2
        center_1_Y = offsetY + sqrt(3)/2 * rad
        center_2_X = offsetX + rad/2
        center_2_Y = offsetY - sqrt(3)/2 * rad
        center_3_X = offsetX - rad
        center_3_Y = offsetY
        
        pgon1 = nsidedpoly(6, 'Center', [center_1_X, center_1_Y], 'Sidelength', rad)
        pgon2 = nsidedpoly(6, 'Center', [center_2_X, center_2_Y], 'Sidelength', rad)
        pgon3 = nsidedpoly(6, 'Center', [center_3_X, center_3_Y], 'Sidelength', rad)    
        
        [user_X1,user_Y1, v_x1, v_y1] = distributeUsers(center_1_X, center_1_Y, rad, 1, 2)
        [user_X2,user_Y2, v_x2, v_y2] = distributeUsers(center_2_X, center_2_Y, rad, 2, 2)    
        [user_X3,user_Y3, v_x3, v_y3] = distributeUsers(center_3_X, center_3_Y, rad, 3, 2)    
        
        users_X = horzcat(users_X, user_X1, user_X2, user_X3)
        users_Y = horzcat(users_Y, user_Y1, user_Y2, user_Y3)
        
        
        plot(v_x1, v_y1,'g.')
        plot(v_x2, v_y2,'g.')
        plot(v_x3, v_y3,'g.')
        
        
        %plot(user_X1, user_Y1,'bl.')
        plot([pgon1, pgon2, pgon3])
    end
    plot(users_X, users_Y,'bl.')
end

function [gx, gy, r_x, r_y] = distributeUsers(xc,yc,rad, numPolygon, towerDist)
    N = 1000; %Number of users
    R = rad; %Radius of Hexagon
    gy=[]
    gx=[]
    v_x = (R * cos((0:6)*pi/3))+xc;
    v_y = (R * sin((0:6)*pi/3))+yc;
    
    c_x = (R-rand(1, 3*N)*2*R)+xc;
    c_y = (R-rand(1, 3*N)*2*R)+yc;

   if numPolygon == 3
       t = 1
       predel1 = 20
       predel2 = 40
   elseif numPolygon == 2
       t = 3
       predel1 = 40
       predel2 = 60
   elseif numPolygon == 1
       t = 5
       predel1 = 0
       predel2 = 20
   end
    
    r_x = (towerDist * cos((predel1:predel2)*pi/30))+v_x(t);
    r_y = (towerDist * sin((predel1:predel2)*pi/30))+v_y(t);
     
    IN = inpolygon(c_x, c_y, v_x, v_y);

    c_x = c_x(IN);
    c_y = c_y(IN);

    idx = randperm(length(c_x));
    c_x = c_x(idx(1:N));
    c_y = c_y(idx(1:N));
    gx = horzcat(gx, c_x)
    gy = horzcat(gy, c_y)
    
    for i = 1:length(gx)
        dist = sqrt((gx(i)-v_x(t))^2 + (gy(i)-v_y(t))^2)
        if dist < towerDist
            gx(i) = xc + 1/(i*randi(10))
            gy(i) = yc + 1/(i*randi(10))
        end
    end
            
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
