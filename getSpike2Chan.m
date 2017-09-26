function [EODR, TEOD, chName] = getSpike2Chan(ADC, ADC_CH)

if ischar(ADC)
    ADC = load(ADC);
end
prefix = getSpike2Prefix(ADC);
sCh = getfield(ADC, sprintf('%s_Ch%d', prefix, ADC_CH));
EODR = sCh.values;
TEOD = (1:numel(EODR)) * sCh.interval;
chName = sCh.title;
EODR(1) = [];
TEOD(1) = [];

if nargout == 0
    figure; plot(TEOD, EODR);
    xlabel('Time (s)');
    ylabel(sprintf('%s (ch%d)', chName, ADC_CH));
    title('Spike2 time-series');
end
end