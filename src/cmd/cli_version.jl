# Lazy version detection — avoids running external `git` at const-initialization
# time, which triggers a Julia 1.12 compiler inference bug.
const version_cache = Ref{String}("")
const hash_cache = Ref{String}("")
const project_toml = abspath(joinpath(@__DIR__, "..", "..", "Project.toml"))

function run_git_or_unknown(args::Cmd)
    Sys.which("git") === nothing && return "unknown"
    buf = IOBuffer()
    proc = run(pipeline(args, stdout=buf, stderr=devnull); wait=false)
    wait(proc)
    success(proc) || return "unknown"
    return strip(String(take!(buf)))
end

"Read `version` from the package Project.toml; empty when unavailable."
function read_project_version()
    isfile(project_toml) || return ""
    for line in eachline(project_toml)
        m = match(r"""^version\s*=\s*"(.+)"\s*$""", line)
        m !== nothing && return m.captures[1]
    end
    return ""
end

function software_version()
    if isempty(version_cache[])
        pkg = read_project_version()
        version_cache[] = if !isempty(pkg)
            pkg
        else
            run_git_or_unknown(`git -C $(dirname(project_toml)) describe --tags --abbrev=0`)
        end
    end
    return version_cache[]
end

function software_git_hash()
    if isempty(hash_cache[])
        hash_cache[] = run_git_or_unknown(`git -C $(dirname(project_toml)) rev-parse HEAD`)
    end
    return hash_cache[]
end

software_version_label() = "$(software_version()) (git $(software_git_hash()))"
