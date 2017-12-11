# **Model Predictive Control**

**Model Predictive Control Project**

1. Clone/fork the project's template files from the project repository. (Note: Please do not submit your project as a pull request against our repo!)
2. Choose your state, input(s), dynamics, constraints and implement MPC.
3. Test your solution on basic examples.
4. Test your solution on the simulator!
5. When the vehicle is able to drive successfully around the track, submit! Remember to include a separate file in .txt, .md, .html, or .pdf format addressing the questions in the rubric.
6. Try to see how fast you get the vehicle to SAFELY go!

[//]: # (Image References)

[image1]: ./25_5_cte.png "25 5 cte"
[image2]: ./25_5_delta.png "25 5 delta"
[image3]: ./10_15_cte.png "10 15 cte"
[image4]: ./10_15_delta.png "10 15 delta"
[image5]: ./1_1_1_1_1_1_1_cte.png "1 1 1 1 1 1 1 cte"
[image6]: ./1_1_1_1_1_1_1_delta.png "1 1 1 1 1 1 1 delta"
[image7]: ./1000_1000_1_1_1_1_1_cte.png "1000 1000 1 1 1 1 1 cte"
[image8]: ./1000_1000_1_1_1_1_1_delta.png "1000 1000 1 1 1 1 1 delta"

## Rubic Points
### Here I will consider the [rubic points](https://review.udacity.com/#!/rubrics/896/view) individually and describe how I addressed each point in my implementation.

---

### Compilation

#### 1. Code must compile without errors with `cmake` and `make`.

My code could be compiled without errors.

To enable debugging, I appended the following line to CMakeLists.txt:

```cmake
set(CMAKE_BUILD_TYPE Debug)
```
 
### Implementation

#### 2. Student describes their model in detail. This includes the state, actuators and update equations.

I defined the vehicle model as the following fomula:

$$x_{t+1} = x_{t} + v_{t} \times \cos(\psi_{t}) \times dt$$
$$y_{t+1} = y_{t} + v_{t} \times \sin(\psi_{t}) \times dt$$
$$\psi_{t+1} = \psi_{t} + (v_{t}/L_{f}) \times \delta_{t} \times dt$$
$$v_{t+1} = v_{t} + a_{t} \times dt$$
$$cte_{t+1} = f(x_{t}) - y_{t} + v_{t} \times \sin(e \psi_{t}) \times dt$$
$$e\psi_{t+1} = \psi_{t} - \psi des_{t} + (v_{t}/L_{f}) \times \delta_{t} \times dt$$

with constraints: 

$$-0.436332 rad < \delta < 0.436332 rad$$
$$and$$
$$-1 < a < 1$$

First, I pass the current state $(x_{1}, y_{1}, \phi_{1}, v_{1}, cte_{1}, e\phi_{1})$ which I will describe later to the model.

Next, the optimization solver is called. The solver uses the initial state, the model constraints and cost function to return a vector of control inputs $(\delta_{i}, a_{i})$ that minimize the cost function. I applied the first control input to the vehicle and repeat the loop.

Below I will explain how I selected the cost function $J$ defined as following:

$$J = \sum_{t=0}^{N-1} (w_{cte}\cdot cte_{t}^{2} + w_{\psi}\cdot e\psi_{t}^{2} + w_{v}\cdot (v_{t} - v_{ref})^{2}) + \sum_{t=0}^{N-2} (w_{\delta}\cdot \delta_{t}^{2} + w_{a}\cdot a_{t}^{2}) + \sum_{t=0}^{t-3} (w_{\Delta\delta}\cdot (\delta_{t+1} - \delta_{t})^2 + w_{\Delta a}\cdot (a_{t+1} - a_{t})^{2})$$

where $w_{X}$ are weight values.

If I set all the weight values to 1, the CTE value became large as shown in the following graph and the behicle was out of the course.

![alt text][image5]
![alt text][image6]

Because I thought that setting larger weight of CTE $(w_{cte})$ and psi $(w_{\psi})$ would reduce the increase in CTE, I set these weights to 1000. However, the results did not improve as shown in the graph below.

![alt text][image7]
![alt text][image8]

Finally I decided to set $w_{\delta}$ to 100 and $w_{\Delta\delta}$ to 1000. As a result, the change in the steering angle became smoothly, and the behicle was able to make one round of the course without meandering.

![alt text][image1]
![alt text][image2]

The following is the final weights I selected.

| weight             | value |
|:------------------:|:-----:|
| $w_{cte}$          | 1     |
| $w_{\psi}$         | 1     |
| $w_{v}$            | 1     |
| $w_{\delta}$       | 100   |
| $w_{a}$            | 1     |
| $w_{\Delta\delta}$ | 1000  |
| $w_{\Delta a}$     | 1     |

#### 3. Student discusses the reasoning behind the chosen `N` and `dt` values.

First I selected the same N and dt values as Lesson 19 (i.e. N = 25 and dt = 0.05). The resulting CTE and steering angle were somewhat wavy as the following graph:

![alt text][image1]
![alt text][image2]

In order to suppress this small amplitude changes, I increased the timestep duration $dt$ to 0.15 and accordingly reduced the number of timestep to 10.

![alt text][image3]
![alt text][image4]

#### 4. A polynomial is fitted to waypoints.

**waypoints**

I converted the positions of the waypoints to the vehicles' coordinate system (main.cpp l.110):

```cpp
vector<double> v_ptsx;
vector<double> v_ptsy;
for (size_t i = 0; i < ptsx.size(); i++) {
    double dx = ptsx[i] - px;
    double dy = ptsy[i] - py;
    v_ptsx.push_back(dx * cos(-psi) - dy * sin(-psi));
    v_ptsy.push_back(dx * sin(-psi) + dy * cos(-psi));
}
```

**the vehicle state**

Because of conversion to the vehicle coordinate system, x, y and phi of the vehicle state are 0. Also, the fitting polynomial, cte and epsi are calcuated with the following code (main.cpp l.119):

```cpp
auto coeffs = polyfit(Eigen::Map<Eigen::VectorXd>(&v_ptsx[0],v_ptsx.size()),
                    Eigen::Map<Eigen::VectorXd>(&v_ptsy[0],v_ptsy.size()), 3);
double cte = polyeval(coeffs, 0);
double epsi = - atan(coeffs[1]);
```

As a point of note, the speed is in mph, I converted it to mps with the following code (main.cpp l.19, l.95):

```cpp
double mph2mps(double v) { return v * 0.44704; }
double v = mph2mps(j[1]["speed"]);
```

#### 5. The student implements Model Predictive Control that handles a 100 millisecond latency.

I predicted the position of the car after 100 milliseconds by the following calculation:

$$dt = 0.1$$
$$x' = x + v \times \cos(\psi) \times dt$$
$$y' = y + v \times \sin(\psi) \times dt$$
$$\psi' = \psi - (v / Lf) \times \delta \times dt$$
$$v' = v + a \times dt$$ 

The following is the actual code (main.cpp l.90):

```cpp
double dt = 0.1;
vector<double> ptsx = j[1]["ptsx"];
vector<double> ptsy = j[1]["ptsy"];
double px = j[1]["x"];
double py = j[1]["y"];
double psi = j[1]["psi"];
double v = mph2mps(j[1]["speed"]);
double steering_angle = j[1]["steering_angle"];
double throttle = j[1]["throttle"];
px += v * cos(psi) * dt;
py += v * sin(psi) * dt;
psi -= (v / MPC::Lf) * steering_angle * dt;
v += throttle * dt;
```

### Simulation

#### 6. The vehicle must successfully drive a lap around the track.

The car did not pop up onto ledges or roll over any surface.