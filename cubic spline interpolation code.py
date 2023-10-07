from matplotlib import pyplot as plt
import numpy as np
import math
import random


def points_to_matrices(points):
    X = np.zeros((4*(len(points)-1), 4*(len(points)-1)))
    Y = np.zeros((4*(len(points)-1), 1))
    for n in range(len(points)-1):
        X[2*n, 4*n] = math.pow(points[n][0], 3)
        X[2*n, 4*n + 1] = math.pow(points[n][0], 2)
        X[2*n, 4*n+2] = points[n][0]
        X[2*n, 4*n+3] = 1
        Y[2*n, 0] = points[n][1]
    for n in range(1, len(points)):
        X[2*n-1, 4*(n-1)] = math.pow(points[n][0], 3)
        X[2*n-1, 4*(n-1)+1] = math.pow(points[n][0], 2)
        X[2*n-1, 4*(n-1)+2] = math.pow(points[n][0], 1)
        X[2*n-1, 4*(n-1)+3] = math.pow(points[n][0], 0)
        Y[2*n-1, 0] = points[n][1]
    for n in range(1, len(points)-1):
        row = 2*(len(points)-1)+n-1
        X[row, 4*n] = -3*math.pow(points[n][0], 2)
        X[row, 4*n+1] = -2*math.pow(points[n][0], 1)
        X[row, 4*n+2] = -1
        X[row, 4*(n-1)] = 3*math.pow(points[n][0], 2)
        X[row, 4*(n-1)+1] = 2*math.pow(points[n][0], 1)
        X[row, 4*(n-1)+2] = 1
    for n in range(1, len(points)-1):
        row = 2*(len(points)-1)+len(points)-1-1-1+n
        X[row, 4*n] = -6*points[n][0]
        X[row, 4*n+1] = -2
        X[row, 4*(n-1)] = 6*points[n][0]
        X[row, 4*(n-1)+1] = 2
    X[-2, 0] = 3*math.pow(points[0][0], 2)
    X[-2, 1] = 2*points[0][0]
    X[-2, 2] = 1
    X[-1, -4] = 3*math.pow(points[-1][0], 2)
    X[-1, -3] = 2*points[-1][0]
    X[-1, -2] = 1

    return X, Y


def coeff_to_func(K, points):
    coeffs = []
    for i in range(len(K)):
        # print(i)
        if (i+1) % 4 == 0:
            # print(K[i-3])
            coeffs.append(K[i-3:i+1])
    x_bounds = []
    for i in range(len(points)-1):
        x_bounds.append((points[i][0], points[i+1][0]))
    return coeffs, x_bounds


def min_max_scaler(coeffs, xbounds, max_acceleration_mss):
    final_accln = 6*coeffs[-1][0]*xbounds[-1][1] + 2*coeffs[-1][1]
    max_accln = final_accln
    x_max_accln = xbounds[-1][1]
    for n_bound in range(len(xbounds)):
        accln = 6*coeffs[n_bound][0]*xbounds[n_bound][0]+2*coeffs[n_bound][1]
        if abs(accln) > abs(max_accln):
            max_accln = accln
            x_max_accln = xbounds[n_bound][0]
    accln_scaler = math.sqrt(max_acceleration_mss/abs(max_accln))

    '''final_vel = 3 * \
        coeffs[-1][0]*math.pow(xbounds[-1][1], 2) + 2 * \
        coeffs[-1][1]*xbounds[-1][1] + coeffs[-1][2]
    max_vel = final_vel
    x_max_vel = xbounds[-1][1]
    for n_bound in range(len(xbounds)):
        vel = 3*coeffs[n_bound][0]*math.pow(xbounds[n_bound][0], 2) + \
            2*coeffs[n_bound][1]*xbounds[n_bound][0] + coeffs[n_bound][2]
        if abs(vel) > abs(max_vel):
            max_vel = vel
            x_max_vel = xbounds[n_bound][0]
    vel_scaler = max_velocity_ms/abs(max_vel)'''
    return accln_scaler  # , vel_scaler, x_max_vel


def get_universal_scaler(pointsAll, limits):
    min_scaler = 2047
    for n_points in range(len(pointsAll)):
        X, Y = points_to_matrices(pointsAll[n_points])
        Xinv = np.linalg.inv(X)
        K = np.matmul(Xinv, Y)
        coeffs, x_bounds = coeff_to_func(K, pointsAll[n_points])
        accln_scaler = min_max_scaler(coeffs, x_bounds, limits[n_points])
        # print(accln_scaler)
        if accln_scaler < min_scaler:
            min_scaler = accln_scaler
    return min_scaler


def cubic_spline_interpolation(points, accln_scaler):
    X, Y = points_to_matrices(points)
    Xinv = np.linalg.inv(X)
    K = np.matmul(Xinv, Y)
    coeffs, x_bounds = coeff_to_func(K, points)
    # print(x_bounds)
    #plt.scatter([point[0] for point in points], [point[1] for point in points])
    #fig, ax = plt.subplots(3, 1, sharex=True)
    #fig.suptitle("Displacement, Velocity, Acceleration")
    why = []
    deewhy = []
    ex = []
    for n_spline in range(len(coeffs)):
        x = np.linspace(x_bounds[n_spline][0]/accln_scaler,
                        x_bounds[n_spline][1]/accln_scaler, 1000)
        # x=x*accln_scaler
        print("{"+str(x_bounds[n_spline][0]/accln_scaler) +
              "<x<"+str(x_bounds[n_spline][1]/accln_scaler)+"}")
        spline = coeffs[n_spline]

        def _func(x):
            #print(f"{spline[0]}*xxx + {spline[1]}*xx + {spline[2]}*x + {spline[3]}")
            print(
                f"{spline[0]*(accln_scaler**3)}*xxx + {(accln_scaler**2)*spline[1]}*xx + {(accln_scaler)*spline[2]}*x + {spline[3]}")
            return (spline[0]*(accln_scaler**3)*[_x*_x*_x for _x in x] + (accln_scaler**2)*spline[1]*[_x*_x for _x in x] + (accln_scaler)*spline[2]*x + spline[3])

        def _funcprime(x):
            #print(f"vel is : 3*{spline[0]}x + 2*{spline[1]}x + {spline[2]})*{vel_scaler}")
            return (3*spline[0]*(accln_scaler**3)*[_x*_x for _x in x] + 2*(accln_scaler**2)*spline[1]*[_x for _x in x] + (accln_scaler)*spline[2])

        def _funcprimeprime(x):
            #print(f"6*{spline[0]}x + 2*{spline[1]}*{accln_scaler}")
            return (6*spline[0]*(accln_scaler**3)*[_x for _x in x] + 2*(accln_scaler**2)*spline[1])
        y = _func(x)

        #ax[0].plot(x, y)
        #ax[0].set_ylabel("Displacement (m)")
        dy = _funcprime(x)
        #ax[1].plot(x, dy)
        #ax[1].set_ylabel("Velocity (m s^-1)")
        ddy = _funcprimeprime(x)
        #ax[2].plot(x, ddy)
        #ax[2].set_ylabel("Acceleration (m s^-2)")
        # ax[2].set_xlabel("Time")
        why.extend(y)
        deewhy.extend(dy)
        ex.extend(x)
        #plt.plot(x, y)
    plt.show()
    return ex, why, deewhy


points = [(0, 0), (1, 2), (2, 4), (3, 9), (4, 25), (5, 26)]
pointsY = (0, 0), (2, 2.236), (4, 4.474), (9,
                                           9.579), (25, 25.618), (26, 39.368)
pointsY = [(point[1], point[0]) for point in pointsY]
pointsX = [(0, 0), (1, 2.236), (2, 4.474),
           (3, 9.579), (4, 25.618), (5, 39.368)]
pointsX = [(point[1], point[0]) for point in pointsX]

scale = get_universal_scaler([pointsX, pointsY], [3.5, 3.5])

xs = cubic_spline_interpolation(pointsX, scale)
ys = cubic_spline_interpolation(pointsY, scale)
#plt.plot(xs, ys)
# plt.show()


'''
import random
def PID(setpoint, kP, kI, kD, current_pos=0, speed=0.1, limit=30):
    prev_error = setpoint
    integral=0
    positions = [current_pos]
    ticker=0
    while(current_pos!=setpoint):
        error = setpoint-current_pos
        integral+=error
        derivative = error-prev_error
        print(error)
        prop = kP*error
        inte = kI*integral
        der = kD*derivative
        #print(der)
        correction = prop+inte+der
        c_speed = speed + correction
        current_pos+=c_speed#+random.gauss(0, .3)
        positions.append(current_pos)
        prev_error=error
        ticker+=1
        if (ticker>limit):
            break
    plt.plot([x for x in range(len(positions))],positions)
    plt.plot([x for x in range(len(positions))],[setpoint for i in range(len(positions))])
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
PID(3, .3, 0.05, 2, 0, 0.0)'''


def pathFollowPID(datas, kP, kI, kD, reaction=0):
    current_vel = 0
    current_pos = 0
    prev_error = 0
    integral = 0
    prior_mu = 0
    prior_sigma = 1
    posterior_mu = 0
    posterior_sigma = 0
    positions = []
    for i in range(len(datas[0])):
        setpoint = datas[2][i]
        cv2 = (current_vel+random.gauss(0, .3))
        current_vel = (current_vel+random.gauss(0, .3))
        posterior_mu = (cv2*0.09 + current_vel*0.09)/(0.18)
        posterior_sigma = (.09*.09/.18)
        #posterior_mu = (posterior_mu*prior_sigma*prior_sigma + prior_mu*posterior_sigma*posterior_sigma)/(posterior_sigma**2 + prior_sigma**2)
        #posterior_sigma = (posterior_sigma**2)*(prior_sigma**2)/(posterior_sigma**2 + prior_sigma**2)
        prior_mu = posterior_mu
        prior_sigma = posterior_sigma
        error = setpoint-posterior_mu
        prop = kP*error
        integral += error
        inte = kI*integral
        derivative = error-prev_error
        dere = kD*derivative
        correction = prop+inte+dere-reaction
        current_vel += correction
        try:
            current_pos += ((datas[0][i+1]-datas[0][i])*current_vel) + \
                random.gauss(0, (datas[0][i+1]-datas[0][i])*.3)
        except:
            current_pos += 0*current_vel
        positions.append(current_pos)
        prev_error = error
    return positions


'''plt.plot(xs[0], pathFollowPID(xs, .1, .4, 0.02, .2))
plt.plot(xs[0], xs[1])

plt.plot(ys[0], pathFollowPID(ys, .1, .4, 0.02, .2))
plt.plot(ys[0], ys[1])'''

robot_x = pathFollowPID(xs, .1, .4, 0.02, .2)
robot_y = pathFollowPID(ys, .1, .4, 0.02, .2)
plt.plot(robot_x, robot_y)
plt.plot(xs[1], ys[1])
