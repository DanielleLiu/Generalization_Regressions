function [eps] = defineShortPositiveSplits(nantype)


names={'Adaptation_{early}',...
    'PosShort_{early}','PosShort2_{early}','PosShort3_{early}',...
   'PosShort4_{early}','PostShort5_{early}',...
   'Adaptation_{late}',...
    'PosShort_{late}','PosShort2_{late}','PosShort3_{late}',...
   'PosShort4_{late}','PostShort5_{late}'};

eps=defineEpochs(names,...
                 {'Adaptation',...
                'Pos Short','Pos Short 2','Pos Short 3',...
                'Pos Short 4','Pos Short 5',...
                'Adaptation',...
                'Pos Short','Pos Short 2','Pos Short 3',...
                'Pos Short 4','Pos Short 5'},...
                [10,...
                10, 10, 10,...
                10 10,...
                -40,...
                -10, -10, -10,...
                -10 -10],...
                [1,...
                1,1,1,...
                1,1,...
                0,...
                0,0,0,...
                0,0],...
                [0,...
                0,0,0,...
                0,0,...
                5,...
                5,5,5,...
                5,5],...
                nantype);