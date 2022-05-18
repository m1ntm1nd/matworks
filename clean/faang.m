hold on
isd = 500;
rad = isd / 3;
numUePerCell = 100;
minDist = 50;
numTiers = 1;
noise_power = 5;
acceptableErrorOfLinearCoordinate = 10;
xhex = [0];
yhex = [0];
basis_towers_coord_X = [0, 500, 250];
basis_towers_coord_Y = [0, 0, 433.0127];

[towers_X, towers_Y, users_X, users_Y] = addTowersRec2(numTiers, isd, rad, minDist, numUePerCell);
towersMatrix = [towers_X(:), towers_Y(:)];
usersMatrix = [users_X(:), users_Y(:), zeros(length(users_X), 3)];
[usersMatrix] = addNearestSt(usersMatrix, towersMatrix);
[usersMatrixPolar] = ToPolar(usersMatrix);
[towersMatrixPolar] = ToPolar(towersMatrix);
towersMatrixPolar(1, 2) = 0;


detectedMatrix = zeros(length(usersMatrix), 2);
len = length(usersMatrixPolar);

for i = 1:len
    [detectedMatrix(i, 1), detectedMatrix(i, 2)] = calcPosRotated(towersMatrixPolar, usersMatrixPolar, usersMatrix, i, noise_power);
    scatter(detectedMatrix(i, 1), detectedMatrix(i, 2));
end

figure(2);

for c = noise_power:10:50
    noise_power = c;
    for i = 1:len
        [detectedMatrix(i, 1), detectedMatrix(i, 2)] = calcPosRotated(towersMatrixPolar, usersMatrixPolar, usersMatrix, i, noise_power);
    end
    [err_users_X, err_users_Y] = calcFangError(usersMatrix, detectedMatrix);
    [vector_T, vector_CDF] = genCDF(err_users_X, err_users_Y, acceptableErrorOfLinearCoordinate);
    plot(vector_T, vector_CDF);
    hold on
end


figure(3)
errorsVector = [];

lenVec = 0:10:50
for c = lenVec
    noise_power = c
    for i = 1:len
        [detectedMatrix(i, 1), detectedMatrix(i, 2)] = calcPosRotated(towersMatrixPolar, usersMatrixPolar, usersMatrix, i, noise_power);
    end
    [err_users_X, err_users_Y] = calcFangError(usersMatrix, detectedMatrix);
    [meanError] = calcMeanError(err_users_X, err_users_Y);
    errorsVector = horzcat(errorsVector, meanError);
end
plot(lenVec, errorsVector);

xlabel('Mean Error (m^2)');
ylabel('Noise power (m or dB)');
function [meanError] = calcMeanError(err_users_X, err_users_Y)
    vectorErrors = [];
    for c = 1:length(err_users_X)
        vectorErrors = horzcat(vectorErrors, sqrt((err_users_X(c))^2 + (err_users_Y(c))^2));
    end
    
    meanError = sum(vectorErrors)/length(vectorErrors);
end

%calculating distance differens (TDoA*c)
%i = userIndex
%towers_coord_X, towers_coord_Y - mas of 3 coords of stations
function [detected_x, detected_y] = calcPosRotated(towersMatrixPolar, usersMatrixPolar, usersMatrix, i, noise_power)
    maxd = 500/2;
    detected_x = 0;
    detected_y = 0;
    theta = pi/6;
    towers_coord_X = [0, 500, 250];
    while theta <= pi
        if (usersMatrixPolar(i, 2) <= theta) || (usersMatrixPolar(i, 2) >= 2*pi - theta)
            if usersMatrixPolar(i, 2) <= pi
                rotate = -(theta-pi/6);
                towers_coord_Y = [0, 0, -433.0127];
            else
                rotate = (theta-pi/6);
                towers_coord_Y = [0, 0, 433.0127];
            end
                %take coords of closest towers in polar system
                %towers_coord_R = towersMatrixPolar(usersMatrix(i,3:5), 1);
                %towers_coord_Phi = towersMatrixPolar(usersMatrix(i,3:5), 2);
                %towers_coord_Phi = towers_coord_Phi + rotate;
                
                %its fake coord rotated by angle
                [user_coord_x, user_coord_y] = ToDecartesDot(usersMatrixPolar(i, 1), usersMatrixPolar(i, 2)+rotate);
                
                %calc towers in Decartes
                %towers_coord_X = zeros(3, 1);
                %towers_coord_Y = zeros(3, 1);
                
                
                %for i = 1:3
                %    [towers_coord_X(i), towers_coord_Y(i)] = ToDecartesDot(towers_coord_R(i), towers_coord_Phi(i));
                %end
                
                distDiff1 = (sqrt((user_coord_x - towers_coord_X(2))^2 + (user_coord_y - towers_coord_Y(2))^2) - sqrt((user_coord_x - towers_coord_X(1))^2 + (user_coord_y - towers_coord_Y(1))^2));
                distDiff2 = (sqrt((user_coord_x - towers_coord_X(3))^2 + (user_coord_y - towers_coord_Y(3))^2) - sqrt((user_coord_x - towers_coord_X(1))^2 + (user_coord_y - towers_coord_Y(1))^2));
                
                %randn - gaussian distr coefficients w/ disp = 1
                distDiff1 = distDiff1 + noise_power * randn(1,1);
                distDiff2 = distDiff2 + noise_power * randn(1,1);
                
                [fake_detected_x, fake_detected_y] = calcPosition(distDiff1, distDiff2, towers_coord_X, towers_coord_Y);
                [fake_r, fake_phi] = ToPolarDot(fake_detected_x, fake_detected_y);
                [detected_x, detected_y] = ToDecartesDot(fake_r, fake_phi-rotate);
                
                
                distFalsy = sqrt(detected_x^2 + detected_y^2);
                if distFalsy > maxd
                    detected_x = usersMatrix(i,1) + noise_power * randn(1,1);
                    detected_y = usersMatrix(i,2) + noise_power * randn(1,1);
                end
                
        end
        theta = theta + pi/6;
    end    
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

function [x, y] = ToDecartesDot(r, phi)
    x = r * cos(phi);
    y = r * sin(phi);
end

function [r, phi] = ToPolarDot(x, y)
    phi = atan(y/x);
    r = sqrt(x^2 + y^2);
    if x>0 && y<=0
        phi = phi + 2 * pi;
    end
    if x<0
        phi = phi + pi;
    end
end


function [usersMatrixPolar] = ToPolar(usersMatrix)
    n = length(usersMatrix);
    usersMatrixPolar = zeros(n, 2);
    for i = 1:n
        x = usersMatrix(i, 1);
        y = usersMatrix(i, 2);
        phi = atan(y/x);
        rad = sqrt(x^2 + y^2);
        if x>0 && y<=0
            phi = phi + 2* pi;
        end
        if x<0
            phi = phi + pi;
        end
        usersMatrixPolar(i, 1) = rad;
        usersMatrixPolar(i, 2) = phi;
    end
end

function [vector_T, vector_CDF] = genCDF(err_users_X, err_users_Y, acceptableErrorOfLinearCoordinate)
    %vector_S = acceptableErrorOfLinearCoordinate ** 2;
    vector_S = [];
    vector_CDF = [0];
    
    for c = 1:length(err_users_X)
        vector_S = horzcat(vector_S, sqrt((err_users_X(c))^2 + (err_users_Y(c))^2));
    end
    max_err = acceptableErrorOfLinearCoordinate ^ 2;
    vector_T = 0:0.1:max_err;
    for max_S_Curr = 0.1:0.1:max_err
        cnt = 0;
        for i = 1:length(vector_S)
            if vector_S(i) <= max_S_Curr
                cnt= cnt+1;
            end
        end
        vector_CDF = horzcat(vector_CDF, cnt/length(err_users_X));       
    end
    
end

function [err_users_X, err_users_Y] = calcFangError(usersMatrix, detectedMatrix)
    err_users_X = [];
    err_users_Y = [];
    
    len = length(usersMatrix);
    for i = 1:len
        if detectedMatrix(i, 1) ~= 0
            err_users_X = horzcat(err_users_X, abs(usersMatrix(i, 1) - detectedMatrix(i, 1)));
            err_users_Y = horzcat(err_users_Y, abs(usersMatrix(i, 2) - detectedMatrix(i, 2)));
        end
    end
end



function [usersMatrix] = addNearestSt(usersMatrix, towersMatrix)
    for i = 1:length(usersMatrix)
        user_X = usersMatrix(i, 1);
        user_Y = usersMatrix(i, 2);   
        [indexes] = findNearestStationsForUser(towersMatrix, user_X, user_Y);
        usersMatrix(i,3:5) = indexes;
    end
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

%util functions
function [indexes] = findNearestStationsForUser(towersMatrix, user_X, user_Y)
    xc = [];
    yc = [];
    rads=[];
    indexes = [];
    for i = 1:length(towersMatrix)
        tower_X = towersMatrix(i, 1);
        tower_Y = towersMatrix(i, 2);   
        rads = horzcat(rads, sqrt((tower_X-user_X)^2 + (tower_Y-user_Y)^2));
    end
    [rads, indexes] = sort(rads);
    indexes = indexes(1:3);
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
