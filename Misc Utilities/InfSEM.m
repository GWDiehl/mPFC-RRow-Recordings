function SEM = InfSEM(input)

usable = isfinite(input);
values = length(find(usable));

SEM = std(input(usable))/sqrt(values);