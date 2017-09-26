function y = semboot(vrY)

nBoot = 1000;
vec = bootci(nBoot, {@(y)mean(y), vrY}, 'type', 'cper');
y = abs(diff(vec)/2);