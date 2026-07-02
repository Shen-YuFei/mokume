use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct StringIdRegistry<I> {
    forward: HashMap<String, I>,
    reverse: Vec<String>,
}

impl<I> Default for StringIdRegistry<I> {
    fn default() -> Self {
        Self {
            forward: HashMap::new(),
            reverse: Vec::new(),
        }
    }
}

impl<I> StringIdRegistry<I>
where
    I: Copy + From<u32> + Into<u32>,
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.reverse.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reverse.is_empty()
    }

    pub fn get(&self, value: &str) -> Option<I> {
        self.forward.get(value).copied()
    }

    pub fn resolve(&self, id: I) -> Option<&str> {
        let raw: u32 = id.into();
        let index = usize::try_from(raw).ok()?;
        self.reverse.get(index).map(|value| value.as_str())
    }

    pub fn get_or_insert(&mut self, value: &str) -> Option<I> {
        if let Some(id) = self.get(value) {
            return Some(id);
        }
        let index = u32::try_from(self.reverse.len()).ok()?;
        let id = I::from(index);
        self.reverse.push(value.to_owned());
        self.forward.insert(value.to_owned(), id);
        Some(id)
    }

    pub fn iter(&self) -> impl Iterator<Item = (I, &str)> + '_ {
        self.reverse
            .iter()
            .enumerate()
            .filter_map(|(index, value)| {
                u32::try_from(index)
                    .ok()
                    .map(|index| (I::from(index), value.as_str()))
            })
    }
}
