
name_list = {'hibiscus'}; % 'bump' 

test_type = 'mesh';
test_type = 'image';


for iname = 1:length(name_list)
    name = name_list{iname};
    for m_porous = [1 1.1 1.2 1.4 1.5 2 4 6]
        disp(['---- ' name ' m_porous=' num2str(m_porous) ' ----']);
        eval(['test_' test_type ';']);
    end
end