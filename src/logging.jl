#this variable will be added
function BotNotify(text::String; add_date = true)
    if add_date
         println("$(now()) $text")
         sendMessage(text = "$(now()) $text")
    else
         println(text)
         sendMessage(text = text) 
    end
    return nothing
end

function BotFigure(fig)
     open(fig, "r") do io
          sendPhoto(photo = io)
      end
end