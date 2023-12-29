abstract type SimCondition end

function integrate1!(atoms::Atoms, condition::SimCondition)
end

function preneighbor!(atoms::Atoms, condition::SimCondition)
end

function integrate2!(atoms::Atoms, condition::SimCondition)
end

function postforce!(atoms::Atoms, condition::SimCondition)
end

function endofstep!(atoms::Atoms, condition::SimCondition)
end
