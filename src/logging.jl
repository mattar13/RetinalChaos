#this variable will be added
function BotNotify(text::String; add_date = true)
    if add_date
         println("[$(now())]: $text")
         sendMessage(text = "[$(now())]: $text")
    else
         println(text)
         sendMessage(text = text) 
    end
    return nothing
end

function BotFigure(fig_loc)
     open(fig_loc, "r") do io
          sendPhoto(photo = io)
      end
end