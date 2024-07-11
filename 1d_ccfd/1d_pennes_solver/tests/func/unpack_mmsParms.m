function [base,startPower,trials,expectedOrder,tolerance,displayPlot,displayData] = unpack_mmsParms(mmsParms)

	if mmsParms == 'default'
		base          = 10; 
		startPower    = 1;
		trials        = 4;
		expectedOrder = 2;
		tolerance     = 0.1;
		displayPlot   = 0;
		displayData   = 0;
	else
		base          = mmsParms{1};
		startPower    = mmsParms{2};
		trials        = mmsParms{3};
		expectedOrder = mmsParms{4};
		tolerance     = mmsParms{5};
		displayPlot   = mmsParms{6};
		displayData   = mmsParms{7};
	end
end
