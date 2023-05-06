clear
clc

%%uniform quantization
% reconstruction from oversampling
t=0:(1/4000)*1:1;% time signal
x = cos(2*pi*2*t);
numlist = [3,4,5,10];
y = zeros(length(numlist),length(t));

for i = 1:length(numlist)
    y(i,:) = double(fi(x,1,2*numlist(i)+1,numlist(i)));
    figure
    plot(t, y(i,:),'r' )
    hold on;
    plot(t, x,'b' )
    title('Quantized signal (8 bits)');
    xlabel('time')
    ylabel('quantized signal')
    hold off
end

error = zeros(1,length(numlist));

for i = 1:length(numlist)
    error(i) = calculate_abs_mean_square_error(x,y(i,:));
end
figure
plot(numlist, error,'r' )

%%non unifrom quantization (mu law)
% reconstruction from oversampling
t=0:(1/4000)*1:1;% time signal
x = cos(2*pi*2*t);
x = ((log(1+255*abs(x)))/(log(1+255))).*sign(x);
numlist = [3,4,5,10];
y = zeros(length(numlist),length(t));

for i = 1:length(numlist)
    y(i,:) = double(fi(x,1,2*numlist(i)+1,numlist(i)));
    figure
    plot(t, y(i,:),'r' )
    hold on;
    plot(t, x,'b' )
    title('Quantized signal (8 bits)');
    xlabel('time')
    ylabel('quantized signal')
    hold off
end

error = zeros(1,length(numlist));

for i = 1:length(numlist)
    error(i) = calculate_abs_mean_square_error(x,y(i,:));
end
figure
plot(numlist, error,'r' )

function abs_mean_square_error = calculate_abs_mean_square_error(vector1, vector2)
% This function calculates the absolute mean square error between two vectors.

% Calculate the difference between the two vectors
difference = vector1 - vector2;

% Square each element of the difference vector
squared_difference = difference .^ 2;

% Calculate the mean of the squared difference vector
mean_squared_difference = mean(squared_difference);

% Calculate the square root of the mean squared difference
abs_mean_square_error = sqrt(mean_squared_difference);
end

