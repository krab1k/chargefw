target_sources(chargefw_core
        PRIVATE
        atom.cpp
        bond.cpp
        conformer.cpp
        molecule.cpp
        molecule_collection.cpp

        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ../../include
        FILES
        ../../include/chargefw/core/atom.h
        ../../include/chargefw/core/bond.h
        ../../include/chargefw/core/conformer.h
        ../../include/chargefw/core/molecule.h
        ../../include/chargefw/core/molecule_collection.h
        ../../include/chargefw/core/position.h
)