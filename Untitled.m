num_sample = 1001; 
nvar = 10;
nobj = 2;
x = zeros(num_sample, nvar);
y = zeros(num_sample, nobj);
result_x = [];
for i = 1 : num_sample
    x(i, 1) = (i - 1) / (num_sample - 1);
    for j = 2 : nvar
        x(i, j) = sin(6 * pi * x(i, 1) + j * pi / nvar);
    end
    y(i,:) = evulate_cf3(x(i,:));
    if constraint(y(i,:)) == 1
        result_x = [result_x; x(i,:)];
    end
end
plot3(result_x(:,1),result_x(:,2),result_x(:,3),'ro');


