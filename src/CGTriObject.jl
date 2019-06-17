# abstract description of a CG elastic object with triangle elements
mutable struct CGTriObject <: ElasticObject
### Attributes
    
### Constructor
    function CGTriObject(
        mesh::Mesh,
        # TODO: Psi::EnergyDensityFunction,
        mass::Float64)
        
    end
end
