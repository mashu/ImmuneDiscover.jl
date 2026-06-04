module UnicodePlotsExt

# Package extension: loaded automatically when UnicodePlots is available. It wires the
# barplot hook so Data.barplot_if_available draws a unicode bar plot; without UnicodePlots
# the hook stays unset and the code falls back to a plain text table.

using UnicodePlots
import immunediscover

function __init__()
    immunediscover.Data.barplot_fn[] = (x, y) -> UnicodePlots.barplot(x, y)
end

end
