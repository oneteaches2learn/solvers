function d = generate_test_data(j)

    % test data
    hand_mean   = [ 32,  27,  25,  22,  20,  18,  18,  18,  17,  16,  15,  16,  15];
    hand_stdDev = [  2,   1,   2,   1,   2,   1,   2,   3,   2,   2,   2,   2,   2];
    fing_mean   = [ 30,  16,  13,  11,  10,  11,  12,  12,  10,  11,  12,  12,  12];
    fing_stdDev = [  3,   4,   4,   2,   2,   2,   2,   2,   4,   2,   2,   2,   3];
    core_mean   = [ 37,  37,  37,  37,  37,  37,  37,  37,  37,  37,  37,  37,  37];
    core_stdDev = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3] / 0.3;

    means     = [hand_mean ; fing_mean ; core_mean];
    stdDevs   = [hand_stdDev ; fing_stdDev ; core_stdDev];

    d = GaussianParameter(hand_mean(j), hand_stdDev(j));

    d = GaussianParameter(means(:,j),stdDevs(:,j));

    d = GaussianParameter(means,stdDevs);

end

