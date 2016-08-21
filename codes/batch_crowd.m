name_list = {'moomoo'};
name_list = { 'aniso0'  'aniso1' 'aniso2' 'aniso3' 'aniso4'};
name_list = {'disk1' 'disk' 'holes' 'tworooms'};
name_list = {'tworooms'};

test_type = 'mesh';
test_type = 'image';


for iname = 1:length(name_list)
    name = name_list{iname};
    for kappa_mult = [1 2 4 6]
        disp(['---- ' name ' kappa_mult=' num2str(kappa_mult) ' ----']);
        eval(['test_' test_type ';']);
    end
end
