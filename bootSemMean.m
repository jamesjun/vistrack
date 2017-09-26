function [sem1] = bootSemMean(vrY)

nBoot = 1000;
vec = bootci(nBoot, {@(y)mean(y), vrY}, 'type', 'cper');
sem1 = diff(vec)/2;