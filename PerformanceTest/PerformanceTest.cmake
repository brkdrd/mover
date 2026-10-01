# Root
set(PERFORMANCE_TEST_ROOT ${PHYSICS_REPO_ROOT}/PerformanceTest)

# Source files
set(PERFORMANCE_TEST_SRC_FILES
	${PERFORMANCE_TEST_ROOT}/PyramidScene.h
	${PERFORMANCE_TEST_ROOT}/PerformanceTest.cpp
	${PERFORMANCE_TEST_ROOT}/PerformanceTest.cmake
	${PERFORMANCE_TEST_ROOT}/PerformanceTestScene.h
	# TODO(4D): these scenes need features that are cut from v1 - Ragdoll, Character, and the
	# Mesh shape (ConvexVsMesh / LargeMesh). Left in the tree, out of the build; restore them
	# with their subsystems.
#	${PERFORMANCE_TEST_ROOT}/RagdollScene.h
#	${PERFORMANCE_TEST_ROOT}/ConvexVsMeshScene.h
#	${PERFORMANCE_TEST_ROOT}/CharacterVirtualScene.h
#	${PERFORMANCE_TEST_ROOT}/LargeMeshScene.h
	${PERFORMANCE_TEST_ROOT}/Layers.h
	${PERFORMANCE_TEST_ROOT}/MaxBodiesScene.h
)

# Group source files
source_group(TREE ${PERFORMANCE_TEST_ROOT} FILES ${PERFORMANCE_TEST_SRC_FILES})
