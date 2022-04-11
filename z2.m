hold on
isd = 500;
rad = isd / 3;
numUePerCell = 10;
minDist = 50;
numTiers = 1;
noise_power = 0.05;
acceptableErrorOfLinearCoordinate = 5;
xhex = [0];
yhex = [0]; 

[towers_X, towers_Y, users_X, users_Y] = addTowersRec2(numTiers, isd, rad, minDist, numUePerCell);
[err_users_X, err_users_Y] = calcFangError(towers_X, towers_Y, users_X, users_Y, noise_power);
[vector_T, vector_CDF] = genCDF(err_users_X, err_users_Y, acceptableErrorOfLinearCoordinate);

figure(2)
plot(vector_T, vector_CDF)

%compare real coordinates with those calculated by the fang algorithm
disp(err_users_X);
%disp(err_users_Y);


function [vector_T, vector_CDF] = genCDF(err_users_X, err_users_Y, acceptableErrorOfLinearCoordinate)
    %vector_S = acceptableErrorOfLinearCoordinate ** 2;
    vector_S = [];
    vector_CDF = [0];
    
    for c = 1:length(err_users_X)
        vector_S = horzcat(vector_S, ((2*err_users_X(c))^2 + (2*err_users_Y(c))^2));
    end
    max_S = acceptableErrorOfLinearCoordinate ^ 2;
    vector_T = 0:max_S;
    for max_S_Curr = 1:max_S
        cnt = 0;
        for i = 1:length(vector_S)
            if vector_S(i) <= max_S_Curr
                cnt= cnt+1;
            end
        end
        vector_CDF = horzcat(vector_CDF, cnt/length(err_users_X));       
    end
    
end

function [err_users_X, err_users_Y] = calcFangError(towers_X, towers_Y, users_X, users_Y, noise_power)
    err_users_X = [];
    err_users_Y = [];

    for c = 1:length(users_X)
        [distDiff1, distDiff2] = calcTDoA(towers_X, towers_Y, users_X(c), users_Y(c), noise_power);
        [xc, yc] = calcPosition(distDiff1, distDiff2, towers_X, towers_Y);
        err_users_X = horzcat(err_users_X, abs(xc - users_X(c)));
        err_users_Y = horzcat(err_users_Y, abs(yc - users_Y(c)));
    end
    
end

%calculating distance differens (TDoA*c)
function [distDiff1, distDiff2] = calcTDoA(towers_coord_X, towers_coord_Y, user_coord_x, user_coord_y, noise_power)
    distDiff1 = abs(sqrt((user_coord_x - towers_coord_X(2))^2 + (user_coord_y - towers_coord_Y(2))^2) - sqrt((user_coord_x - towers_coord_X(1))^2 + (user_coord_y - towers_coord_Y(1))^2));
    distDiff2 = abs(sqrt((user_coord_x - towers_coord_X(3))^2 + (user_coord_y - towers_coord_Y(3))^2) - sqrt((user_coord_x - towers_coord_X(1))^2 + (user_coord_y - towers_coord_Y(1))^2));
    
    %add white gaussian noise *2 (on transmitter and receiver)
    distDiff1 = awgn(distDiff1, 2/noise_power);
    distDiff2 = awgn(distDiff2, 2/noise_power);
end

function [xc, yc] = calcPosition(distDiff1, distDiff2, towers_coord_X, towers_coord_Y)
    g = (distDiff2*(towers_coord_X(2)/distDiff1) - towers_coord_X(3))/towers_coord_Y(3);
    h = ((((towers_coord_X(3))^2 + (towers_coord_Y(3))^2)) - distDiff2^2 + distDiff2 * distDiff1*(1-(towers_coord_X(2)/distDiff1)^2)) / (2*towers_coord_Y(3));

    d = -(1-(towers_coord_X(2)/distDiff1)^2 + g^2);
    e = towers_coord_X(2)*(1-(towers_coord_X(2)/distDiff1)^2)-(2*g*h);
    f = ((distDiff1^2) / 4)*(1-(towers_coord_X(2)/distDiff1)^2)^2 - h^2;

    xc = (-e-sqrt(e^2 - 4*d*f))/(2*d);
    yc = g * xc + h;
end




function [towers_coord_X,towers_coord_Y, users_X, users_Y] = addTowersRec2(numTiers, isd, rad, minDist, numUePerCell)
    towers_coord_X = [0];
    towers_coord_Y = [0];    
    x=[1 1/2 -1/2 -1 -1/2 1/2] * isd;
    y=[0 sqrt(3)/2 sqrt(3)/2 0 -sqrt(3)/2 -sqrt(3)/2] * isd;
    for newTire = 1:numTiers
        for c = 1:length(towers_coord_X)
            offsetX = towers_coord_X(c);
            offsetY = towers_coord_Y(c);
            new_towers_group_coords_X = x - offsetX;
            new_towers_group_coords_Y = y - offsetY;
            [new_towers_X, new_towers_Y] = check_unique_coords(towers_coord_X, towers_coord_Y, new_towers_group_coords_X, new_towers_group_coords_Y);
            towers_coord_X = horzcat(towers_coord_X, new_towers_X);
            towers_coord_Y = horzcat(towers_coord_Y, new_towers_Y);
        end
    end

    plot(towers_coord_X,towers_coord_Y, 'Or');

    users_X = [];
    users_Y = [];

    DISTR_USERS_FLAG = 0;
    for c = 1:length(towers_coord_X)
        offsetX = towers_coord_X(c);
        offsetY = towers_coord_Y(c);
        center_1_X = offsetX + rad/2;
        center_1_Y = offsetY + sqrt(3)/2 * rad;
        center_2_X = offsetX + rad/2;
        center_2_Y = offsetY - sqrt(3)/2 * rad;
        center_3_X = offsetX - rad;
        center_3_Y = offsetY;

        pgon1 = nsidedpoly(6, 'Center', [center_1_X, center_1_Y], 'Sidelength', rad);
        pgon2 = nsidedpoly(6, 'Center', [center_2_X, center_2_Y], 'Sidelength', rad);
        pgon3 = nsidedpoly(6, 'Center', [center_3_X, center_3_Y], 'Sidelength', rad);    
        
        if DISTR_USERS_FLAG == 0
            [user_X1,user_Y1, v_x1, v_y1] = distributeUsers(center_1_X, center_1_Y, rad, 1, minDist, numUePerCell);
            [user_X2,user_Y2, v_x2, v_y2] = distributeUsers(center_2_X, center_2_Y, rad, 2, minDist, numUePerCell);    
            [user_X3,user_Y3, v_x3, v_y3] = distributeUsers(center_3_X, center_3_Y, rad, 3, minDist, numUePerCell);    
       
        
        
            users_X = horzcat(users_X, user_X1, user_X2, user_X3);
            users_Y = horzcat(users_Y, user_Y1, user_Y2, user_Y3);
            DISTR_USERS_FLAG = 1;
        end

        plot(v_x1, v_y1,'g.');
        plot(v_x2, v_y2,'g.');
        plot(v_x3, v_y3,'g.');


        %plot(user_X1, user_Y1,'bl.');
        plot([pgon1, pgon2, pgon3]);
    end
    plot(users_X, users_Y,'bl.');
end

function [gx, gy, r_x, r_y] = distributeUsers(xc,yc,rad, numPolygon, towerDist, numUePerCell)
    if towerDist > rad * 2
       error('Error. towerDist can not be larger than the diameter.');
    end

    gy=[];
    gx=[];

    v_x = (rad * cos((0:6)*pi/3))+xc;
    v_y = (rad * sin((0:6)*pi/3))+yc;

    if numPolygon == 3
        t = 1;
        predel1 = 20;
        predel2 = 40;
    elseif numPolygon == 2
        t = 3;
        predel1 = 40;
        predel2 = 60;
    elseif numPolygon == 1
        t = 5;
        predel1 = 0;
        predel2 = 20;
    end

    r_x = (towerDist * cos((predel1:predel2)*pi/30))+v_x(t);
    r_y = (towerDist * sin((predel1:predel2)*pi/30))+v_y(t);

    while length(gx) < numUePerCell
        c_x = (rad - rand(1, 3*numUePerCell)*2*rad)+xc;
        c_y = (rad - rand(1, 3*numUePerCell)*2*rad)+yc;

        IN = inpolygon(c_x, c_y, v_x, v_y);

        c_x = c_x(IN);
        c_y = c_y(IN);

        idx = randperm(length(c_x));
        c_x = c_x(idx(1:numUePerCell));
        c_y = c_y(idx(1:numUePerCell));

        for i = 1:length(c_x)
            dist = sqrt((c_x(i)-v_x(t))^2 + (c_y(i)-v_y(t))^2);
            if (dist > towerDist) && (length(gx) < numUePerCell)
                gx = horzcat(gx, c_x(i));
                gy = horzcat(gy, c_y(i));
            end
        end
    end      
end

function [mas_X, mas_Y ] = check_unique_coords(towers_coord_X, towers_coord_Y, new_towers_group_coords_X, new_towers_group_coords_Y)
    mas_X = [];
    mas_Y = [];
    for i = 1:length(new_towers_group_coords_X)
        flag = 0;
        curr_check_dot_X = new_towers_group_coords_X(i);
        curr_check_dot_Y = new_towers_group_coords_Y(i);
        for j = 1:length(towers_coord_X)
            if ((curr_check_dot_X == towers_coord_X(j)) && (curr_check_dot_Y == towers_coord_Y(j)))
                flag = 1;
            end
        end

        if (flag == 0)
            mas_X = horzcat(mas_X, curr_check_dot_X);
            mas_Y = horzcat(mas_Y, curr_check_dot_Y);
        end
    end
end
