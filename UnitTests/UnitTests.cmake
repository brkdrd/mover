# Root
set(UNIT_TESTS_ROOT ${PHYSICS_REPO_ROOT}/UnitTests)

# Source files
set(UNIT_TESTS_SRC_FILES
#	${UNIT_TESTS_ROOT}/Compute/ComputeTests.cpp	# TODO(4D): still 3D (Mat44 SetColumn3/Multiply3x3), rewrite for 4D
	${UNIT_TESTS_ROOT}/Core/ArrayTest.cpp
	${UNIT_TESTS_ROOT}/Core/BinaryHeapTest.cpp
	${UNIT_TESTS_ROOT}/Core/FPFlushDenormalsTest.cpp
	${UNIT_TESTS_ROOT}/Core/HashCombineTest.cpp
	${UNIT_TESTS_ROOT}/Core/InsertionSortTest.cpp
	${UNIT_TESTS_ROOT}/Core/JobSystemTest.cpp
	${UNIT_TESTS_ROOT}/Core/LinearCurveTest.cpp
	${UNIT_TESTS_ROOT}/Core/PreciseMathTest.cpp
	${UNIT_TESTS_ROOT}/Core/ScopeExitTest.cpp
	${UNIT_TESTS_ROOT}/Core/STLLocalAllocatorTest.cpp
	${UNIT_TESTS_ROOT}/Core/StringToolsTest.cpp
	${UNIT_TESTS_ROOT}/Core/QuickSortTest.cpp
	${UNIT_TESTS_ROOT}/Core/UnorderedSetTest.cpp
	${UNIT_TESTS_ROOT}/Core/UnorderedMapTest.cpp
	${UNIT_TESTS_ROOT}/doctest.h
	${UNIT_TESTS_ROOT}/Geometry/AABoxTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/ClipPolytopeTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/ClosestPointTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/ConvexHullBuilderTest.cpp
	${UNIT_TESTS_ROOT}/Geometry/EllipseTest.cpp
	${UNIT_TESTS_ROOT}/Geometry/EPATests.cpp
	${UNIT_TESTS_ROOT}/Geometry/GJKTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/PlaneTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/RayAABoxTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/RaySphereTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/RayTetrahedronTests.cpp
	${UNIT_TESTS_ROOT}/Geometry/RayTriangleTests.cpp
	${UNIT_TESTS_ROOT}/Layers.h
	${UNIT_TESTS_ROOT}/LoggingBodyActivationListener.h
	${UNIT_TESTS_ROOT}/LoggingContactListener.h
	${UNIT_TESTS_ROOT}/Math/BVec16Tests.cpp
	${UNIT_TESTS_ROOT}/Math/BivecTests.cpp
	${UNIT_TESTS_ROOT}/Math/Cl4Reference.h
	${UNIT_TESTS_ROOT}/Math/DVec4Tests.cpp
	# TODO(4D): these tests are still written against the 3D API (Vec3/Quat, the cut
	# Mesh/HeightField/Cylinder/Tapered*/Triangle shapes, joint constraints, SoftBody, Character,
	# Vehicle) and need rewriting for 4D before they can be re-registered. Re-enable them one
	# by one as each area is ported - do NOT delete, they encode a lot of expected behaviour.
#	${UNIT_TESTS_ROOT}/Math/EigenValueSymmetricTests.cpp
	${UNIT_TESTS_ROOT}/Math/HalfFloatTests.cpp
	${UNIT_TESTS_ROOT}/Math/Mat44DeterminantTests.cpp
#	${UNIT_TESTS_ROOT}/Math/Mat44Tests.cpp	# TODO(4D): still 3D (Mat44 SetColumn3/Multiply3x3), rewrite for 4D
	${UNIT_TESTS_ROOT}/Math/MathTests.cpp
	${UNIT_TESTS_ROOT}/Math/MatrixTests.cpp
#	${UNIT_TESTS_ROOT}/Math/QuatTests.cpp
	${UNIT_TESTS_ROOT}/Math/RMat44Tests.cpp
	${UNIT_TESTS_ROOT}/Math/RotorTests.cpp
	${UNIT_TESTS_ROOT}/Math/UVec4Tests.cpp
	${UNIT_TESTS_ROOT}/Math/TrigonometryTests.cpp
	${UNIT_TESTS_ROOT}/Math/Vec3Tests.cpp
	${UNIT_TESTS_ROOT}/Math/Vec4Tests.cpp
	${UNIT_TESTS_ROOT}/Math/Vec8Tests.cpp
	${UNIT_TESTS_ROOT}/Math/VectorTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/ActiveEdgesTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/BroadPhaseTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/CastShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/CharacterVirtualTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/CollideShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/CollidePointTests.cpp
	${UNIT_TESTS_ROOT}/Physics/CollisionGroupTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/ContactListenerTests.cpp
	${UNIT_TESTS_ROOT}/Physics/ConstraintPartTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/ConvexVsTrianglesTest.cpp
#	${UNIT_TESTS_ROOT}/Physics/DistanceConstraintTests.cpp
	${UNIT_TESTS_ROOT}/Physics/DynamicsConventionsTest.cpp
#	${UNIT_TESTS_ROOT}/Physics/EstimateCollisionResponseTest.cpp
#	${UNIT_TESTS_ROOT}/Physics/HeightFieldShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/HingeConstraintTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/MotionQualityLinearCastTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/MutableCompoundShapeTests.cpp
	${UNIT_TESTS_ROOT}/Physics/ObjectLayerPairFilterTableTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/ObjectLayerPairFilterMaskTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/OffsetCenterOfMassShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/PathConstraintTests.cpp	# TODO(4D): PathConstraint is cut for v1
#	${UNIT_TESTS_ROOT}/Physics/PhysicsDeterminismTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/PhysicsStepListenerTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/PhysicsTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/RayShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/SensorTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/ShapeFilterTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/ShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/SixDOFConstraintTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/SliderConstraintTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/SoftBodyTests.cpp
	${UNIT_TESTS_ROOT}/Physics/SubShapeIDTest.cpp
#	${UNIT_TESTS_ROOT}/Physics/TaperedCylinderShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/TransformedShapeTests.cpp
#	${UNIT_TESTS_ROOT}/Physics/WheeledVehicleTests.cpp
#	${UNIT_TESTS_ROOT}/PhysicsTestContext.cpp
	${UNIT_TESTS_ROOT}/PhysicsTestContext.h
	${UNIT_TESTS_ROOT}/UnitTestFramework.cpp
	${UNIT_TESTS_ROOT}/UnitTestFramework.h
	${UNIT_TESTS_ROOT}/UnitTests.cmake
)

if (ENABLE_OBJECT_STREAM)
	set(UNIT_TESTS_SRC_FILES
		${UNIT_TESTS_SRC_FILES}
#		${UNIT_TESTS_ROOT}/ObjectStream/ObjectStreamTest.cpp
		${UNIT_TESTS_ROOT}/ObjectStream/RotorBivecStreamTest.cpp
	)
endif()

# Assets used by the unit tests
set(UNIT_TESTS_ASSETS
	${JOLT_PHYSICS_SPV_SHADERS}
	${JOLT_PHYSICS_METAL_LIB}
)

# Group source files
source_group(TREE ${UNIT_TESTS_ROOT} FILES ${UNIT_TESTS_SRC_FILES})
