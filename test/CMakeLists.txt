add_executable(Action_test "Action_test.cpp")
target_link_libraries( Action_test Actions)

add_executable( Integrator_test "Integrator_test.cpp" )
target_link_libraries( Integrator_test HMC_Support )

add_executable(Metropolis_test "Metropolis_test.cpp")
target_link_libraries( Metropolis_test UpdateAlgorithms)

add_executable(HMC_test "HMC_test.cpp")
target_link_libraries( HMC_test UpdateAlgorithms)