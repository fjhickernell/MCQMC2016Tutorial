%ut_workout_integral_g unit test for workout_integral_g
classdef ut_workout_integral_g < matlab.unittest.TestCase

  methods(Test)
             
    function testWorkout_integral_g(testCase)
      nrep=10000; nmax=1e7; abstol=1e-8;
      [succnowarn, succwarn, pfin] = workout_integral_g(nrep,nmax,abstol);
      succrates = succnowarn + succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,pfin);
      testCase.verifyGreaterThan(succrates(3), succrates(2));
      testCase.verifyGreaterThan(succrates(2), succrates(1));
    end
    
  end
end