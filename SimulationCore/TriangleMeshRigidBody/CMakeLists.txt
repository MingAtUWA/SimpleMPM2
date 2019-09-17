target_sources(
    SimulationCore
    PRIVATE
        # Rigid body
        Geometry.h Geometry.cpp
        HashTable.h HashTableEdge.h
        TriangleMesh.h TriangleMesh.cpp TriangleMesh_BgGrid.cpp
        RigidBody.h RigidBody.cpp
    )