function SEM = SEM(input)

usable = isfinite(input);
values = length(find(usable));

SEM = std(input(usable))/sqrt(values);