export force!, ZeroAtomic

struct ZeroAtomic <: AtomicPotential end

function force!(atoms::Atoms, params::ZeroAtomic) end
